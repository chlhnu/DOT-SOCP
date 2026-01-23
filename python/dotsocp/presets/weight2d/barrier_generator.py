from abc import ABC, abstractmethod
from pathlib import Path
from typing import Callable

import matplotlib.image as mpimg
import numpy as np
from scipy.interpolate import RegularGridInterpolator

from dotsocp.presets.resources import _resolve_resources_path

__all__ = [
    "get_barrier_fn",
    "get_predefined_barriers",
    "make_barrier_from_img",
    "make_barrier_from_formula",
]


# =========================
# Helpers
# =========================

def _prepare_points(points: np.ndarray) -> tuple[np.ndarray, tuple[int, ...], bool]:
    pts = np.asarray(points, dtype=np.float64)
    squeezed = False
    if pts.ndim == 1:
        if pts.size != 2:
            raise ValueError("points vector must have exactly two entries")
        pts = pts.reshape(1, 2)
        squeezed = True
    elif pts.shape[-1] != 2:
        raise ValueError("points array must have shape (..., 2)")

    value_shape = pts.shape[:-1]
    pts_flat = pts.reshape(-1, 2)
    return pts_flat, value_shape, squeezed


def _circle_pillar_formula(x: np.ndarray, y: np.ndarray) -> np.ndarray:
    left_pillar = (x >= 0.2) & (x <= 0.25) & (y >= 0.4) & (y <= 1.0)
    right_pillar = (x >= 0.75) & (x <= 0.8) & (y >= 0.0) & (y <= 0.6)
    circle = ((x - 0.5) ** 2 + (y - 0.5) ** 2) <= 0.15 ** 2
    mask = left_pillar | right_pillar | circle
    return mask.astype(np.float64)


def _love_heart_formula(x: np.ndarray, y: np.ndarray) -> np.ndarray:
    center1 = np.array([0.7, 0.3], dtype=np.float64)
    center2 = np.array([0.345, 0.625], dtype=np.float64)
    r1 = 0.09
    r2 = 0.09

    sigma1 = r1 / 3.0
    sigma2 = r2 / 3.0

    dist1_sq = (x - center1[0]) ** 2 + (y - center1[1]) ** 2
    dist2_sq = (x - center2[0]) ** 2 + (y - center2[1]) ** 2

    rho0 = np.exp(-dist1_sq / (2.0 * sigma1 ** 2))
    rho0 = np.where(dist1_sq <= r1 ** 2, rho0, 0.0)

    rho1 = np.exp(-dist2_sq / (2.0 * sigma2 ** 2))
    rho1 = np.where(dist2_sq <= r2 ** 2, rho1, 0.0)

    return rho0 + rho1

# =========================
# Factories
# =========================


class BarrierFactory(ABC):
    @abstractmethod
    def create(self) -> Callable[[np.ndarray], np.ndarray]:
        raise NotImplementedError


class BarrierFactoryFromImg(BarrierFactory):
    def __init__(self, image_path: str | Path, method: str = "nearest"):
        self.image_path = Path(_resolve_resources_path(str(image_path)))
        self.method = method

    def create(self) -> Callable[[np.ndarray], np.ndarray]:
        barrier = self._load_barrier()
        interpolator = self._build_interpolator(barrier)
        return self._wrap_interpolator(interpolator)

    def _load_barrier(self) -> np.ndarray:
        barrier = mpimg.imread(self.image_path)
        if barrier.ndim == 3:
            barrier = barrier[..., 0]
        barrier = np.asarray(barrier, dtype=np.float64)

        max_val = float(np.max(barrier))
        if max_val <= 0.0:
            raise ValueError("Barrier image must have positive pixel values")

        return (max_val - barrier) * (255.0 / max_val)

    def _build_interpolator(self, barrier: np.ndarray) -> RegularGridInterpolator:
        nx, ny = barrier.shape
        grid_x = np.linspace(0.0, 1.0, nx)
        grid_y = np.linspace(0.0, 1.0, ny)

        return RegularGridInterpolator(
            (grid_x, grid_y),
            barrier,
            method=self.method,
            bounds_error=False,
            fill_value=None,
        )

    @staticmethod
    def _wrap_interpolator(interpolator: RegularGridInterpolator) -> Callable[[np.ndarray], np.ndarray]:
        def barrier_fn(points: np.ndarray) -> np.ndarray:
            pts_flat, value_shape, squeezed = _prepare_points(points)
            values = interpolator(pts_flat)
            if squeezed:
                return float(values[0])
            return values.reshape(value_shape)

        return barrier_fn


class BarrierFactoryFromFormula(BarrierFactory):
    def __init__(self, formula: Callable[[np.ndarray, np.ndarray], np.ndarray]):
        self.formula = formula

    def create(self) -> Callable[[np.ndarray], np.ndarray]:
        def barrier_fn(points: np.ndarray) -> np.ndarray:
            pts_flat, value_shape, squeezed = _prepare_points(points)
            pts = pts_flat.reshape(value_shape + (2,))
            x = pts[..., 0]
            y = pts[..., 1]
            values = np.asarray(self.formula(x, y), dtype=np.float64)
            try:
                values = np.broadcast_to(values, value_shape)
            except ValueError as exc:
                raise ValueError("Barrier formula must broadcast to the input grid shape") from exc
            if squeezed:
                return float(values.reshape(-1)[0])
            return values.reshape(value_shape)

        return barrier_fn


# =========================
# Public API
# =========================

def make_barrier_from_img(image_name: str, method: str = "nearest") -> Callable[[np.ndarray], np.ndarray]:
    factory = BarrierFactoryFromImg(image_name, method=method)
    return factory.create()

def make_barrier_from_formula(formula: Callable[[np.ndarray, np.ndarray], np.ndarray]) -> Callable[[np.ndarray], np.ndarray]:
    factory = BarrierFactoryFromFormula(formula)
    return factory.create()


# Predefined barriers
_BARRIER_REGISTRY: dict[str, Callable[[], Callable[[np.ndarray], np.ndarray]]] = {
    "love-heart": lambda: make_barrier_from_formula(_love_heart_formula), # Example 5.8: barrier in the shape of a love heart.
    "maze": lambda: make_barrier_from_img("maze.png"), # Example 5.9: barrier in the shape of a maze.
    "maze14": lambda: make_barrier_from_img("maze14.png"), # barrier used in [Optimal Transport with Proximal Splitting. SIAM Journal on Imaging Sciences, 2014.]
    "circle-pillar": lambda: make_barrier_from_formula(_circle_pillar_formula),
}


def get_barrier_fn(name: str) -> Callable[[np.ndarray], np.ndarray]:
    func = _BARRIER_REGISTRY.get(name.lower().replace('_', '-'))
    if func is None:
        raise ValueError(f"Unknown barrier preset: {name}")
    return func()

def get_predefined_barriers() -> list[str]:
    return list(_BARRIER_REGISTRY.keys())

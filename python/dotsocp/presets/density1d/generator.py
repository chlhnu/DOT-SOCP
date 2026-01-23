from typing import Callable
import numpy as np


__all__ = [
    "get_predefined_examples", # List of predefined examples
    "get_example", # Predefined examples
]

# =========================
# Helpers
# =========================

def _grid(nx: int) -> np.ndarray:
    return np.linspace(0.0, 1.0, nx, dtype=np.float64).reshape(nx, 1)


def _normalize_1d(rho: np.ndarray, nx: int, lower_bound: float) -> np.ndarray:
    total = float(np.sum(rho))
    if total <= 0.0:
        raise ValueError("Density must have positive mass for normalization")
    scale = nx / total
    return (scale * rho + lower_bound) / (1.0 + lower_bound)

# =========================
# Example generators
# =========================

def gene_example_gaussian(nx: int) -> tuple[np.ndarray, np.ndarray]:
    mu1 = 0.3
    mu2 = 0.7
    sigma1 = 0.01
    sigma2 = sigma1 / 4.0

    def _normal(x, mean_value, sigma_inv):
        return np.sqrt(sigma_inv) / (2.0 * np.pi) * np.exp(-0.5 * sigma_inv * (x - mean_value) ** 2)

    x = _grid(nx)
    rho0 = _normal(x, mu1, 1.0 / sigma1)
    rho1 = _normal(x, mu2, 1.0 / sigma2)
    return rho0, rho1


def gene_example_box(nx: int) -> tuple[np.ndarray, np.ndarray]:
    x = _grid(nx)
    lbox = (0.1, 0.5)
    rbox = (0.85, 0.95)
    rho0 = ((x >= lbox[0]) & (x <= lbox[1])).astype(np.float64)
    rho1 = ((x >= rbox[0]) & (x <= rbox[1])).astype(np.float64)
    return rho0, rho1

# =========================
# Public API
# =========================

_EXAMPLE_REGISTRY: dict[str, Callable[[int], tuple[np.ndarray, np.ndarray]]] = {
    "gaussian": gene_example_gaussian,
    "box": gene_example_box,
}


def get_example(problem_name: str, nx: int, lower_bound: float = 0.0):
    func = _EXAMPLE_REGISTRY.get(problem_name.lower().replace('_', '-'))
    if func is None:
        # print(f"Warning: Problem '{problem_name}' not implemented, using example 'gaussian'.")
        # func = gene_example_gaussian
        raise ValueError(f"Problem '{problem_name}' not supported yet.")

    rho0, rho1 = func(nx)
    rho0 = _normalize_1d(rho0, nx, lower_bound)
    rho1 = _normalize_1d(rho1, nx, lower_bound)

    return np.asfortranarray(rho0), np.asfortranarray(rho1)


def get_predefined_examples():
    return list(_EXAMPLE_REGISTRY.keys())

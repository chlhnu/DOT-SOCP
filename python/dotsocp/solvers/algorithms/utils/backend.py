from __future__ import annotations

import importlib
import importlib.util
import logging
from dataclasses import dataclass, field
from functools import lru_cache
from types import ModuleType
from typing import Literal

import numpy as np

from dotsocp.configs.icons import LOG_ICON

BackendMode = Literal["auto", "python", "cpp"]

_SQRT_HALF = np.sqrt(2.0) / 2.0
_LOAD_ERRORS: dict[str, BaseException] = {}
_PY_BACKEND_WARNED: set[int] = set()


def _normalize_backend_mode(backend_mode: str) -> BackendMode:
    mode = backend_mode.lower()
    if mode not in {"auto", "python", "cpp"}:
        raise ValueError("backend_mode must be one of: 'auto', 'python', 'cpp'.")
    return mode  # type: ignore[return-value]


def _cpp_module_name(dim: int) -> str:
    if dim == 1:
        return f"{__package__}.cpp_ext_1d"
    if dim == 2:
        return f"{__package__}.cpp_ext_2d"
    raise ValueError(f"Unsupported dimension for C++ backend: {dim}")


@lru_cache(maxsize=2)
def _load_cpp_module(full_name: str) -> ModuleType | None:
    try:
        spec = importlib.util.find_spec(full_name)
    except (ImportError, ValueError) as exc:
        _LOAD_ERRORS[full_name] = exc
        return None

    if spec is None:
        return None

    try:
        return importlib.import_module(full_name)
    except (ImportError, OSError) as exc:
        _LOAD_ERRORS[full_name] = exc
        logging.debug("C++ backend %s is unavailable: %s", full_name, exc)
        return None


def compiled_backend_available(dim: int) -> bool:
    return _load_cpp_module(_cpp_module_name(dim)) is not None


def _warn_python_backend(dim: int) -> None:
    if dim in _PY_BACKEND_WARNED:
        return
    logging.warning("%s  cpp backend for %sD is unavailable; using py backend.\n", LOG_ICON['warn'], dim)
    _PY_BACKEND_WARNED.add(dim)


def _resolve_cpp_module(*, dim: int, backend_mode: str) -> ModuleType | None:
    mode = _normalize_backend_mode(backend_mode)
    if mode == "python":
        return None

    full_name = _cpp_module_name(dim)
    module = _load_cpp_module(full_name)
    if module is not None:
        return module

    if mode == "cpp":
        exc = _LOAD_ERRORS.get(full_name)
        msg = (
            f"C++ backend for {dim}D is not available. Rebuild the extension for "
            "the current Python/platform, or use backend_mode='auto'/'python'."
        )
        raise RuntimeError(msg) from exc

    return None


def _as_f_vector(value: np.ndarray) -> np.ndarray:
    return np.asarray(value, dtype=np.float64).reshape(-1, order="F")


def _require_shape(name: str, value: np.ndarray, shape: tuple[int, int]) -> None:
    if value.shape != shape:
        raise ValueError(f"Expected {name} of shape {shape}, got {value.shape}.")


def _require_size(name: str, value: np.ndarray, size: int) -> None:
    if value.size != size:
        raise ValueError(f"Expected {name} of size {size}, got {value.size}.")


def _write_vector(dst: np.ndarray, src: np.ndarray) -> None:
    dst[...] = np.asarray(src, dtype=np.float64).reshape(dst.shape, order="F")


def _oper_bfd_1d_py(
    z: np.ndarray,
    q: np.ndarray,
    nt: int,
    nx: int,
    scale_bf: float,
    scale_d: float,
) -> None:
    len_z = (nt - 1) * nx
    len_q = len_z + nt * (nx - 1)
    _require_shape("z", z, (len_z, 6))

    q_vec = _as_f_vector(q)
    _require_size("q", q_vec, len_q)

    out = np.zeros((len_z, 6), dtype=np.float64, order="F")
    a = q_vec[:len_z]
    out[:, 0] = scale_d - scale_bf * a
    out[:, 5] = scale_d + scale_bf * a

    b = (scale_bf * _SQRT_HALF) * q_vec[len_z:].reshape((nx - 1, nt), order="F")
    zero = np.zeros((1, nt - 1), dtype=np.float64, order="F")
    out[:, 1] = np.concatenate((zero, b[:, :-1]), axis=0).reshape(-1, order="F")
    out[:, 2] = np.concatenate((b[:, :-1], zero), axis=0).reshape(-1, order="F")
    out[:, 3] = np.concatenate((zero, b[:, 1:]), axis=0).reshape(-1, order="F")
    out[:, 4] = np.concatenate((b[:, 1:], zero), axis=0).reshape(-1, order="F")
    z[...] = out


def _oper_bfd_conj_1d_py(
    q: np.ndarray,
    z: np.ndarray,
    nt: int,
    nx: int,
    scale_bf: float,
) -> None:
    len_z = (nt - 1) * nx
    len_q = len_z + nt * (nx - 1)
    _require_size("q", q, len_q)
    _require_shape("z", z, (len_z, 6))

    a = scale_bf * (z[:, 5] - z[:, 0])
    tmp = z[:, 1:5].reshape((nx, nt - 1, 4), order="F")
    zero = np.zeros((nx - 1, 1), dtype=np.float64, order="F")
    b = (scale_bf * _SQRT_HALF) * (
        np.concatenate((tmp[1:, :, 0], zero), axis=1)
        + np.concatenate((tmp[:-1, :, 1], zero), axis=1)
        + np.concatenate((zero, tmp[1:, :, 2]), axis=1)
        + np.concatenate((zero, tmp[:-1, :, 3]), axis=1)
    )
    _write_vector(q, np.concatenate((a, b.reshape(-1, order="F"))))


def _oper_bfd_2d_py(
    z: np.ndarray,
    q: np.ndarray,
    nt: int,
    nx: int,
    ny: int,
    scale_bf: float,
    scale_d: float,
) -> None:
    len_z = (nt - 1) * nx * ny
    len_b = nt * (nx - 1) * ny
    len_c = nt * nx * (ny - 1)
    len_q = len_z + len_b + len_c
    _require_shape("z", z, (len_z, 10))

    q_vec = _as_f_vector(q)
    _require_size("q", q_vec, len_q)

    out = np.zeros((len_z, 10), dtype=np.float64, order="F")
    a = q_vec[:len_z]
    out[:, 0] = scale_d - scale_bf * a
    out[:, 9] = scale_d + scale_bf * a

    b_start = len_z
    b_end = b_start + len_b
    b = (scale_bf * _SQRT_HALF) * q_vec[b_start:b_end].reshape((ny, nx - 1, nt), order="F")
    zero_x = np.zeros((ny, 1, nt - 1), dtype=np.float64, order="F")
    out[:, 1] = np.concatenate((zero_x, b[:, :, :-1]), axis=1).reshape(-1, order="F")
    out[:, 2] = np.concatenate((b[:, :, :-1], zero_x), axis=1).reshape(-1, order="F")
    out[:, 3] = np.concatenate((zero_x, b[:, :, 1:]), axis=1).reshape(-1, order="F")
    out[:, 4] = np.concatenate((b[:, :, 1:], zero_x), axis=1).reshape(-1, order="F")

    c = (scale_bf * _SQRT_HALF) * q_vec[b_end:].reshape((ny - 1, nx, nt), order="F")
    zero_y = np.zeros((1, nx, nt - 1), dtype=np.float64, order="F")
    out[:, 5] = np.concatenate((zero_y, c[:, :, :-1]), axis=0).reshape(-1, order="F")
    out[:, 6] = np.concatenate((c[:, :, :-1], zero_y), axis=0).reshape(-1, order="F")
    out[:, 7] = np.concatenate((zero_y, c[:, :, 1:]), axis=0).reshape(-1, order="F")
    out[:, 8] = np.concatenate((c[:, :, 1:], zero_y), axis=0).reshape(-1, order="F")
    z[...] = out


def _oper_bfd_conj_2d_py(
    q: np.ndarray,
    z: np.ndarray,
    nt: int,
    nx: int,
    ny: int,
    scale_bf: float,
) -> None:
    len_z = (nt - 1) * nx * ny
    len_b = nt * (nx - 1) * ny
    len_c = nt * nx * (ny - 1)
    len_q = len_z + len_b + len_c
    _require_size("q", q, len_q)
    _require_shape("z", z, (len_z, 10))

    a = scale_bf * (z[:, 9] - z[:, 0])

    tmp_x = z[:, 1:5].reshape((ny, nx, nt - 1, 4), order="F")
    zero_x = np.zeros((ny, nx - 1, 1), dtype=np.float64, order="F")
    b = (scale_bf * _SQRT_HALF) * (
        np.concatenate((tmp_x[:, 1:, :, 0], zero_x), axis=2)
        + np.concatenate((tmp_x[:, :-1, :, 1], zero_x), axis=2)
        + np.concatenate((zero_x, tmp_x[:, 1:, :, 2]), axis=2)
        + np.concatenate((zero_x, tmp_x[:, :-1, :, 3]), axis=2)
    )

    tmp_y = z[:, 5:9].reshape((ny, nx, nt - 1, 4), order="F")
    zero_y = np.zeros((ny - 1, nx, 1), dtype=np.float64, order="F")
    c = (scale_bf * _SQRT_HALF) * (
        np.concatenate((tmp_y[1:, :, :, 0], zero_y), axis=2)
        + np.concatenate((tmp_y[:-1, :, :, 1], zero_y), axis=2)
        + np.concatenate((zero_y, tmp_y[1:, :, :, 2]), axis=2)
        + np.concatenate((zero_y, tmp_y[:-1, :, :, 3]), axis=2)
    )

    _write_vector(
        q,
        np.concatenate(
            (
                a,
                b.reshape(-1, order="F"),
                c.reshape(-1, order="F"),
            )
        ),
    )


def _proj_soc_py(z: np.ndarray, zz: np.ndarray) -> None:
    if z.shape != zz.shape:
        raise ValueError(f"Expected z and zz to share the same shape, got {z.shape} and {zz.shape}.")
    if zz.ndim != 2 or zz.shape[1] < 2:
        raise ValueError("proj_soc expects 2-D arrays with at least two columns.")

    zz_arr = np.asarray(zz, dtype=np.float64)
    t = zz_arr[:, 0]
    x = zz_arr[:, 1:]
    x_norm = np.linalg.norm(x, axis=1)
    out = np.zeros_like(zz_arr, dtype=np.float64, order="F")

    mask_inside = x_norm <= t
    out[mask_inside] = zz_arr[mask_inside]

    mask_mid = (~mask_inside) & (x_norm > -t) & (x_norm > 0.0)
    if np.any(mask_mid):
        factor = 0.5 * (1.0 + t[mask_mid] / x_norm[mask_mid])
        out[mask_mid, 0] = factor * x_norm[mask_mid]
        out[mask_mid, 1:] = factor[:, None] * x[mask_mid]

    z[...] = out


@dataclass(slots=True)
class SocpBackend:
    dim: int
    nt: int
    nx: int
    ny: int | None = None
    backend_mode: str = "auto"
    _cpp_module: ModuleType | None = field(init=False, default=None, repr=False)

    def __post_init__(self) -> None:
        self.backend_mode = _normalize_backend_mode(self.backend_mode)
        if self.dim == 2 and self.ny is None:
            raise ValueError("ny must be provided for a 2D backend.")
        self._cpp_module = _resolve_cpp_module(dim=self.dim, backend_mode=self.backend_mode)
        if self._cpp_module is None and self.backend_mode == "auto":
            _warn_python_backend(self.dim)

    @property
    def using_compiled_backend(self) -> bool:
        return self._cpp_module is not None

    def oper_bfd(
        self,
        z: np.ndarray,
        q: np.ndarray,
        scale_bf: float = 1.0,
        scale_d: float = 1.0,
    ) -> None:
        if self._cpp_module is not None:
            if self.dim == 1:
                self._cpp_module.oper_bfd(z, q, self.nt, self.nx, float(scale_bf), float(scale_d))
            else:
                assert self.ny is not None
                self._cpp_module.oper_bfd(z, q, self.nt, self.nx, self.ny, float(scale_bf), float(scale_d))
            return

        if self.dim == 1:
            _oper_bfd_1d_py(z, q, self.nt, self.nx, float(scale_bf), float(scale_d))
            return

        assert self.ny is not None
        _oper_bfd_2d_py(z, q, self.nt, self.nx, self.ny, float(scale_bf), float(scale_d))

    def oper_bfd_conj(
        self,
        q: np.ndarray,
        z: np.ndarray,
        scale_bf: float = 1.0,
    ) -> None:
        if self._cpp_module is not None:
            if self.dim == 1:
                self._cpp_module.oper_bfd_conj(q, z, self.nt, self.nx, float(scale_bf))
            else:
                assert self.ny is not None
                self._cpp_module.oper_bfd_conj(q, z, self.nt, self.nx, self.ny, float(scale_bf))
            return

        if self.dim == 1:
            _oper_bfd_conj_1d_py(q, z, self.nt, self.nx, float(scale_bf))
            return

        assert self.ny is not None
        _oper_bfd_conj_2d_py(q, z, self.nt, self.nx, self.ny, float(scale_bf))

    def proj_soc(self, z: np.ndarray, zz: np.ndarray) -> None:
        if self._cpp_module is not None:
            self._cpp_module.proj_soc(z, zz)
            return
        _proj_soc_py(z, zz)


def make_socp_backend(
    *,
    dim: int,
    nt: int,
    nx: int,
    ny: int | None = None,
    backend_mode: str = "auto",
) -> SocpBackend:
    return SocpBackend(dim=dim, nt=nt, nx=nx, ny=ny, backend_mode=backend_mode)


__all__ = ["SocpBackend", "compiled_backend_available", "make_socp_backend"]

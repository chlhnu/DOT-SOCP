import numpy as np


def _evaluate_barrier(barrier_fn, x_grid, y_grid):
    try:
        values = barrier_fn(x_grid, y_grid)
    except TypeError:
        points = np.stack((x_grid, y_grid), axis=-1)
        values = barrier_fn(points)
    return np.asarray(values, dtype=np.float64)


def _normalize_density(rho: np.ndarray) -> np.ndarray:
    total = float(np.sum(rho))
    if total <= 0.0:
        raise ValueError("Density must have positive total mass after applying barrier constraints")
    scale = rho.size / total
    return rho * scale


def ensure_barrier_validity(rho0: np.ndarray, rho1: np.ndarray, barrier_fn):
    """Project densities so they vanish on the barrier and renormalize mass."""
    rho0 = np.asarray(rho0, dtype=np.float64).copy()
    rho1 = np.asarray(rho1, dtype=np.float64).copy()

    ny, nx = rho0.shape
    if rho1.shape != (ny, nx):
        raise ValueError("rho0 and rho1 must share the same shape")

    x = np.linspace(0.0, 1.0, nx)
    y = np.linspace(0.0, 1.0, ny)
    xx, yy = np.meshgrid(x, y, indexing="xy")

    barrier_grid = _evaluate_barrier(barrier_fn, xx.T, yy.T)
    barrier_grid = np.asarray(barrier_grid, dtype=np.float64)
    if barrier_grid.shape == (nx, ny):
        barrier_grid = barrier_grid.T
    elif barrier_grid.shape != (ny, nx):
        raise ValueError("Barrier grid must match density shape")

    threshold = float(np.mean(barrier_grid))
    barrier_mask = barrier_grid > threshold

    rho0[barrier_mask] = 0.0
    rho1[barrier_mask] = 0.0

    rho0 = _normalize_density(rho0)
    rho1 = _normalize_density(rho1)

    return rho0, rho1


def get_weight_from_barrier(barrier_fn, nt, nx, ny, barrier_weight=1e6):
    """Generate weights from a barrier function."""
    def _eval_barrier(x_vals, y_vals):
        try:
            return np.asarray(barrier_fn(x_vals, y_vals), dtype=np.float64)
        except TypeError:
            pts = np.stack((x_vals, y_vals), axis=-1)
            return np.asarray(barrier_fn(pts), dtype=np.float64)

    hx = 1.0 / (nx - 1)
    hy = 1.0 / (ny - 1)

    x_stag = np.linspace(0.5 * hx, 1 - 0.5 * hx, nx - 1)
    x_cent = np.linspace(0.0, 1.0, nx)
    y_stag = np.linspace(0.5 * hy, 1 - 0.5 * hy, ny - 1)
    y_cent = np.linspace(0.0, 1.0, ny)

    xx, yy = np.meshgrid(x_stag, y_cent, indexing="xy")
    mask = _eval_barrier(xx.T, yy.T) > 0
    weight_x = np.ones((ny, nx - 1))
    weight_x[mask.T] = barrier_weight

    xx, yy = np.meshgrid(x_cent, y_stag, indexing="xy")
    mask = _eval_barrier(xx.T, yy.T) > 0
    weight_y = np.ones((ny - 1, nx))
    weight_y[mask.T] = barrier_weight

    weight_t = np.ones((ny, nx, nt - 1))
    weight_x = np.repeat(weight_x[:, :, None], nt, axis=2)
    weight_y = np.repeat(weight_y[:, :, None], nt, axis=2)

    return np.concatenate([
        weight_t.ravel(order="F"),
        weight_x.ravel(order="F"),
        weight_y.ravel(order="F"),
    ])


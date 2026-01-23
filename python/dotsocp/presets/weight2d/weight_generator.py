import numpy as np


__all__ = [
    "get_weight",
    "get_predefined_weights",
]

# =========================
# Helpers
# =========================

def gene_weight_from_fn(nt, nx, ny, weight_fn):
    """Helper to generate weights from an arbitrary spatial function."""
    hx = 1.0 / (nx - 1)
    hy = 1.0 / (ny - 1)
    x_stag = np.linspace(0.5 * hx, 1 - 0.5 * hx, nx - 1)
    x_cent = np.linspace(0.0, 1.0, nx)
    y_stag = np.linspace(0.5 * hy, 1 - 0.5 * hy, ny - 1)
    y_cent = np.linspace(0.0, 1.0, ny)

    xx, yy = np.meshgrid(x_stag, y_cent)
    weight_x = weight_fn(xx, yy)
    weight_x *= ny * (nx - 1) / np.sum(weight_x)

    xx, yy = np.meshgrid(x_cent, y_stag)
    weight_y = weight_fn(xx, yy)
    weight_y *= ny * (nx - 1) / np.sum(weight_y)

    weight_t = np.ones((ny, nx, nt - 1))
    weight_x = np.tile(weight_x[:, :, None], (1, 1, nt))
    weight_y = np.tile(weight_y[:, :, None], (1, 1, nt))

    return np.concatenate(
        (
            weight_t.ravel(order='F'),
            weight_x.ravel(order='F'),
            weight_y.ravel(order='F'),
        )
    )

# =========================
# Example generators
# =========================


def gene_weight_circular(nt, nx, ny):
    a, b = 0.5, 0.5
    circular = lambda xx, yy: np.sqrt((xx - a) ** 2 + (yy - b) ** 2)
    return gene_weight_from_fn(nt, nx, ny, circular)


def gene_weight_circular_inv(nt, nx, ny):
    a, b = 0.5, 0.5
    circular_inv = lambda xx, yy: 0.1 / (0.1 + np.sqrt((xx - a) ** 2 + (yy - b) ** 2))
    return gene_weight_from_fn(nt, nx, ny, circular_inv)

# =========================
# Public API
# =========================

_WEIGHT_REGISTRY: dict[str, callable] = {
    "circular": gene_weight_circular,
    "circular-inv": gene_weight_circular_inv,
}

def get_weight(weight_name, nt, nx, ny):
    """Get a weight array from the registry."""
    func = _WEIGHT_REGISTRY.get(weight_name.lower().replace('_', '-'))
    if func is None:
        # print(f"Warning: Weight '{weight_name}' not implemented, using weight 'circular'.")
        # func = gene_weight_circular
        raise ValueError(f"Problem '{weight_name}' not supported yet.")

    return func(nt, nx, ny)

def get_predefined_weights():
    """Get all predefined weights."""
    return list(_WEIGHT_REGISTRY.keys())

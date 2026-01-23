import numpy as np
from scipy.fft import dctn, idctn


# =========================
# Opeartors for q variable
# =========================

def oper_q_2d(ny, nx, nt, D=1.0, E=1.0, weight=None):
    if weight is None:
        length = ny * nx * (nt - 1) + ny * (nx - 1) * nt + (ny - 1) * nx * nt
        weight_square = np.broadcast_to(1.0, (length,))
    else:
        weight_square = np.pow(weight, 2).flatten(order='F')

    _scaled_const = (E / D) ** 2
    val_main = 2 * _scaled_const
    val_edge = _scaled_const

    # shape: (ny, nx, nt-1)
    a = np.full((ny, nx, nt - 1), val_main, dtype=np.float64, order='F')

    # shape: (ny, nx-1, nt)
    b = np.full((ny, nx - 1, nt), val_main, dtype=np.float64, order='F')

    b[:, :, 0] = val_edge
    b[:, :, -1] = val_edge

    # shape: (ny-1, nx, nt)
    c = np.full((ny - 1, nx, nt), val_main, dtype=np.float64, order='F')
    c[:, :, 0] = val_edge
    c[:, :, -1] = val_edge

    return np.concatenate([
        a.flatten(order='F'),
        b.flatten(order='F'),
        c.flatten(order='F')
    ]) + weight_square


def oper_q_1d(nx, nt, D=1.0, E=1.0, weight=None):
    if weight is None:
        length = nx * (nt - 1) + (nx - 1) * nt
        weight_square = np.broadcast_to(1.0, (length,))
    else:
        weight_square = np.pow(weight, 2).flatten(order='F')

    _scaled_const = (E / D) ** 2
    val_main = 2.0 * _scaled_const
    val_edge = _scaled_const
    
    a = np.full((nx, nt - 1), val_main, dtype=np.float64, order='F')
    b = np.full((nx - 1, nt), val_main, dtype=np.float64, order='F')
    b[:, 0] = val_edge
    b[:, -1] = val_edge

    return np.concatenate([
        a.flatten(order='F'),
        b.flatten(order='F')
    ]) + weight_square


# ========================
# FFT-based Poisson solver
# =======================

def create_fft_kernel_2d(nt, nx, epsilon=0.0):
    t_idx = np.arange(nt)
    CT = (2 * (nt - 1) ** 2) * (1 - np.cos(np.pi * t_idx / nt))

    x_idx = np.arange(nx)
    CX = (2 * (nx - 1) ** 2) * (1 - np.cos(np.pi * x_idx / nx))

    kernel = CX.reshape(-1, 1) + CT.reshape(1, -1)

    if epsilon > 0:
        kernel += epsilon

    kernel[kernel == 0] = 1.0

    return np.asfortranarray(kernel)


def create_fft_kernel_3d(nt, nx, ny, epsilon=0.0):
    t_idx = np.arange(nt)
    CT = (2 * (nt - 1) ** 2) * (1 - np.cos(np.pi * t_idx / nt))

    x_idx = np.arange(nx)
    CX = (2 * (nx - 1) ** 2) * (1 - np.cos(np.pi * x_idx / nx))

    y_idx = np.arange(ny)
    CY = (2 * (ny - 1) ** 2) * (1 - np.cos(np.pi * y_idx / ny))

    kernel = (CY.reshape(-1, 1, 1) +
              CX.reshape(1, -1, 1) +
              CT.reshape(1, 1, -1))

    if epsilon > 0:
        kernel += epsilon

    kernel[kernel == 0] = 1.0

    return np.asfortranarray(kernel)


def solve_poisson(kernel, rhs):
    assert kernel.size == rhs.size, f"kernel.size ({kernel.size}) != rhs.size ({rhs.size})."

    rhs_nd = rhs.reshape(kernel.shape, order='F')
    rhs_hat = dctn(rhs_nd, type=2, norm='ortho')
    u_hat = rhs_hat / kernel
    u_nd = idctn(u_hat, type=2, norm='ortho')
    u_nd = np.asarray(u_nd)

    return u_nd.reshape(-1, 1, order='F')

from typing import Dict, Optional
from dataclasses import dataclass
import numpy as np
import scipy.sparse as sp


@dataclass
class VariableState():
    qInd: Dict[str, int]
    phi: np.ndarray
    z: np.ndarray
    beta: np.ndarray
    q: np.ndarray
    alpha: np.ndarray
    cScale: Optional[float] = None
    dScale: Optional[float] = None
    D: Optional[float] = None
    E: Optional[float] = None
    E2: Optional[float] = None


@dataclass
class ModelState():
    dim: int
    nx: int
    nt: int
    grad: sp.csr_matrix
    c: np.ndarray
    rho0: np.ndarray
    rho1: np.ndarray
    ht: float
    hx: float
    ny: Optional[int] = None
    hy: Optional[float] = None
    normc: Optional[float] = None
    normd: Optional[float] = None
    weight: Optional[np.ndarray] = None


# =========================
# 1D Initialization
# =========================

def _gene_dt_1d(nt, nx, ht):
    e = np.ones(nt)
    data = np.vstack([-e / ht, e / ht])
    Dt = sp.spdiags(data, [0, 1], nt - 1, nt)
    Ix = sp.eye(nx)
    gradt = sp.kron(Dt, Ix, format='csr')
    return gradt


def _gene_dx_1d(nt, nx, hx):
    It = sp.eye(nt)
    e = np.ones(nx)
    data = np.vstack([-e / hx, e / hx])
    Dx = sp.spdiags(data, [0, 1], nx - 1, nx)
    gradx = sp.kron(It, Dx, format='csr')
    return gradx


def initialize_1d(rho0, rho1, nt, weight=None):
    """Initialize variables and model for 1D DOT-SOCP problem."""

    # Initialize model state
    rho0_arr = np.asarray(rho0, dtype=np.float64).reshape(-1, order='F')
    rho1_arr = np.asarray(rho1, dtype=np.float64).reshape(-1, order='F')

    if rho0_arr.shape[0] != rho1_arr.shape[0]:
        raise ValueError("rho0 and rho1 must have the same length for 1D initialization")

    nx = rho0_arr.shape[0]
    n = nx * nt

    rho0_vec = rho0_arr.reshape(-1, 1, order='F')
    rho1_vec = rho1_arr.reshape(-1, 1, order='F')

    ht = 1.0 / (nt - 1)
    hx = 1.0 / (nx - 1)

    gradt = _gene_dt_1d(nt, nx, ht)
    gradx = _gene_dx_1d(nt, nx, hx)
    grad = sp.vstack([gradt, gradx]).tocsr()

    c = np.zeros((n, 1), dtype=np.float64, order='F')
    c[:nx] = -rho0_vec / ht
    c[-nx:] = rho1_vec / ht

    model = ModelState(
        dim=1,
        nx=nx, nt=nt,
        rho0=rho0, rho1=rho1,
        ht=ht, hx=hx,
        grad=grad,
        c=c,
        weight=weight,
    )

    # Initialize variable state
    len_z = (nt - 1) * nx
    qInd = {'bx': len_z}

    x_lin = np.linspace(0.0, 1.0, nx)
    phi_1d = 0.5 * x_lin ** 2
    phi = np.tile(phi_1d.reshape(-1, 1, order='F'), (1, nt))
    phi = phi.flatten(order='F').reshape(-1, 1)

    z = np.zeros((len_z, 6), dtype=np.float64, order='F')
    beta = np.zeros((len_z, 6), dtype=np.float64, order='F')

    num_edges = grad.shape[0]
    q = np.zeros((num_edges, 1), dtype=np.float64, order='F')
    alpha = np.zeros((num_edges, 1), dtype=np.float64, order='F')

    var = VariableState(
        qInd=qInd,
        phi=phi,
        z=z,
        beta=beta,
        q=q,
        alpha=alpha
    )

    return var, model


# =========================
# 2D Initialization
# =========================

def _gene_dt_2d(nt, nx, ny, ht):
    e = np.ones(nt)
    data = np.vstack([-e / ht, e / ht])
    diags = [0, 1]
    Dt = sp.spdiags(data, diags, nt - 1, nt)
    Ixy = sp.eye(nx * ny)
    gradt = sp.kron(Dt, Ixy, format='csr')
    return gradt


def _gene_dx_2d(nt, nx, ny, hx):
    It = sp.eye(nt)
    Iy = sp.eye(ny)
    e = np.ones(nx)
    data = np.vstack([-e / hx, e / hx])
    Dx = sp.spdiags(data, [0, 1], nx - 1, nx)
    gradx = sp.kron(sp.kron(It, Dx), Iy, format='csr')
    return gradx


def _gene_dy_2d(nt, nx, ny, hy):
    Itx = sp.eye(nt * nx)
    e = np.ones(ny)
    data = np.vstack([-e / hy, e / hy])
    Dy = sp.spdiags(data, [0, 1], ny - 1, ny)
    grady = sp.kron(Itx, Dy, format='csr')
    return grady


def initialize_2d(rho0, rho1, nt, weight=None):
    """Initialize variables and model for 2D DOT-SOCP problem."""
    ny, nx = rho0.shape
    n = nx * ny * nt

    # Initialize model state
    rho0_vec = rho0.flatten(order='F').reshape(-1, 1)
    rho1_vec = rho1.flatten(order='F').reshape(-1, 1)

    ht = 1.0 / (nt - 1)
    hx = 1.0 / (nx - 1)
    hy = 1.0 / (ny - 1)

    gradt = _gene_dt_2d(nt, nx, ny, ht)
    gradx = _gene_dx_2d(nt, nx, ny, hx)
    grady = _gene_dy_2d(nt, nx, ny, hy)

    grad = sp.vstack([gradt, gradx, grady], format='csr')

    c = np.zeros((n, 1), order='F')
    len_space = nx * ny
    c[:len_space] = -rho0_vec / ht
    c[-len_space:] = rho1_vec / ht

    model = ModelState(
        dim=2,
        nx=nx, ny=ny, nt=nt,
        rho0=rho0, rho1=rho1,
        ht=ht, hx=hx, hy=hy,
        grad=grad,
        c=c,
        weight=weight,
    )

    # Initialize variable state
    len_z = (nt - 1) * nx * ny
    idx_bx = len_z
    len_bx = nt * (nx - 1) * ny
    idx_by = idx_bx + len_bx

    qInd = {
        'bx': idx_bx,
        'by': idx_by
    }

    x_lin = np.linspace(0, 1, nx)
    y_lin = np.linspace(0, 1, ny)
    xx, yy = np.meshgrid(x_lin, y_lin, indexing='xy')

    phi_2d = 0.5 * (xx ** 2 + yy ** 2)
    phi_3d = np.repeat(phi_2d[:, :, np.newaxis], nt, axis=2)
    phi = phi_3d.flatten(order='F').reshape(-1, 1)
    z = np.zeros((len_z, 10), order='F')
    beta = np.zeros((len_z, 10), order='F')
    num_edges = grad.shape[0]
    q = np.zeros((num_edges, 1), order='F')
    alpha = np.zeros((num_edges, 1), order='F')

    var = VariableState(
        qInd=qInd,
        phi=phi,
        z=z,
        beta=beta,
        q=q,
        alpha=alpha
    )

    return var, model
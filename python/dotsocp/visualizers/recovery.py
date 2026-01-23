'''Helper functions for recovering the density, movement, and q variables from the solver output.
'''

import numpy as np
from dotsocp.solvers import SolverOutput


def recover_density(result: SolverOutput):
    """Recover density from alpha variables."""
    dim = result.model.dim
    if dim == 1:
        return _recover_density_1d(result)
    if dim == 2:
        return _recover_density_2d(result)
    raise ValueError(f"Unsupported dimension: {dim}")

def recover_movement(result: SolverOutput):
    dim = result.model.dim
    if dim == 1:
        return _recover_movement_1d(result)
    if dim == 2:
        return _recover_movement_2d(result)
    raise ValueError(f"Unsupported dimension: {dim}")

def recover_q(result: SolverOutput):
    dim = result.model.dim
    if dim == 1:
        return _recover_q_1d(result)
    if dim == 2:
        return _recover_q_2d(result)
    raise ValueError(f"Unsupported dimension: {dim}")


def _recover_density_1d(result: SolverOutput):
    alpha = result.var.alpha
    bx_idx = result.var.qInd["bx"]
    nx = result.model.nx
    nt = result.model.nt
    rho0 = result.model.rho0.reshape((nx, 1), order="F")
    rho1 = result.model.rho1.reshape((nx, 1), order="F")

    alpha_rho_vec = alpha[:bx_idx].reshape((nx, nt - 1), order="F")
    transport_rho = np.concatenate((rho0, .5 * alpha_rho_vec[:, :-1] + .5 * alpha_rho_vec[:, 1:], rho1), axis=-1)

    return transport_rho


def _recover_density_2d(result: SolverOutput):
    alpha = result.var.alpha
    bx_idx = result.var.qInd["bx"]
    nx, ny, nt = result.model.nx, result.model.ny, result.model.nt
    rho0 = result.model.rho0.reshape((ny, nx, 1), order="F")
    rho1 = result.model.rho1.reshape((ny, nx, 1), order="F")

    alpha_rho_vec = alpha[:bx_idx].reshape((ny, nx, nt - 1), order="F")
    transport_rho = np.concatenate((rho0, .5 * alpha_rho_vec[:, :, :-1] + .5 * alpha_rho_vec[:, :, 1:], rho1), axis=-1)

    return transport_rho


def _recover_movement_1d(result: SolverOutput):
    """Recover movement (Ex) from alpha on a staggered grid (1D only)."""
    if result.model.dim != 1:
        raise ValueError("recover_movement is only supported for 1D/2D problems.")

    nt = result.model.nt
    nx = result.model.nx
    qInd = result.var.qInd
    alpha = np.asarray(result.var.alpha).reshape(-1, order="F")

    Ex = alpha[qInd["bx"] :].reshape((nx - 1, nt), order="F")
    Ex[:, (0, -1)] *= 2.0
    Ex = np.concatenate(
        (
            np.zeros((1, nt), dtype=Ex.dtype),
            0.5 * (Ex[:-1, :] + Ex[1:, :]),
            np.zeros((1, nt), dtype=Ex.dtype),
        ),
        axis=0,
    )
    return (Ex)


def _recover_movement_2d(result: SolverOutput):
    """
    Recover movement (Ex, Ey) from alpha on a staggered grid (2D only).

    Returns
    -------
    Ex  : ndarray, shape (ny, nx+1, nt-1)
    Ey  : ndarray, shape (ny+1, nx, nt-1)
    """
    if result.model.dim != 2:
        raise ValueError("recover_RhoE is only supported for 2D problems.")

    nt = result.model.nt
    nx = result.model.nx
    ny = result.model.ny

    qInd = result.var.qInd
    alpha = np.asarray(result.var.alpha).reshape(-1, order="F")

    # Ex
    Ex = alpha[qInd["bx"] : qInd["by"]].reshape((ny, nx - 1, nt), order="F")
    Ex[:, :, (0, -1)] *= 2.0
    Ex = np.concatenate(
        (
            np.zeros((ny, 1, nt), dtype=Ex.dtype),
            0.5 * (Ex[:, :-1, :] + Ex[:, 1:, :]),
            np.zeros((ny, 1, nt), dtype=Ex.dtype),
        ),
        axis=1,
    )

    # Ey
    Ey = alpha[qInd["by"] :].reshape((ny - 1, nx, nt), order="F")
    Ey[:, :, (0, -1)] *= 2.0
    Ey = np.concatenate(
        (
            np.zeros((1, nx, nt), dtype=Ey.dtype),
            0.5 * (Ey[:-1, :, :] + Ey[1:, :, :]),
            np.zeros((1, nx, nt), dtype=Ey.dtype),
        ),
        axis=0,
    )

    return (Ex, Ey)


def _recover_q_1d(result: SolverOutput):
    """
    Recover (q0, bx) from q on a time-staggered grid (1D only).

    Returns
    -------
    q0 : ndarray, shape (nx, nt-1)
    bx : ndarray, shape (nx, nt-1)
    """
    if result.model.dim != 1:
        raise ValueError("recover_q is only supported for 1D/2D problems.")

    nt = result.model.nt
    nx = result.model.nx
    qInd = result.var.qInd
    q = np.asarray(result.var.q).reshape(-1, order="F")

    q0 = q[: qInd["bx"]].reshape((nx, nt - 1), order="F")

    bx_core = q[qInd["bx"] :].reshape((nx - 1, nt), order="F")
    bx = np.concatenate(
        (
            np.zeros((1, nt), dtype=bx_core.dtype),
            0.5 * (bx_core[:-1, :] + bx_core[1:, :]),
            np.zeros((1, nt), dtype=bx_core.dtype),
        ),
        axis=0,
    )
    bx = 0.5 * (bx[:, :-1] + bx[:, 1:])

    return (q0, bx)


def _recover_q_2d(result: SolverOutput):
    """
    Recover (q0, bx, by) from q on a time-staggered grid (2D only).

    Returns
    -------
    q0 : ndarray, shape (ny, nx, nt-1)
    bx : ndarray, shape (ny, nx+1, nt-1)
    by : ndarray, shape (ny+1, nx, nt-1)
    """
    if result.model.dim != 2:
        raise ValueError("recover_q is only supported for 2D problems.")

    nt = result.model.nt
    nx = result.model.nx
    ny = result.model.ny

    qInd = result.var.qInd
    q = np.asarray(result.var.q).reshape(-1, order="F")

    # q0 component
    q0 = q[: qInd["bx"]].reshape((ny, nx, nt - 1), order="F")

    # bx component
    bx_core = q[qInd["bx"] : qInd["by"]].reshape((ny, nx - 1, nt), order="F")
    bx = np.concatenate(
        (
            np.zeros((ny, 1, nt), dtype=bx_core.dtype),
            0.5 * (bx_core[:, :-1, :] + bx_core[:, 1:, :]),
            np.zeros((ny, 1, nt), dtype=bx_core.dtype),
        ),
        axis=1,
    )
    bx = 0.5 * (bx[:, :, :-1] + bx[:, :, 1:])

    # by component
    by_core = q[qInd["by"] :].reshape((ny - 1, nx, nt), order="F")
    by = np.concatenate(
        (
            np.zeros((1, nx, nt), dtype=by_core.dtype),
            0.5 * (by_core[:-1, :, :] + by_core[1:, :, :]),
            np.zeros((1, nx, nt), dtype=by_core.dtype),
        ),
        axis=0,
    )
    by = 0.5 * (by[:, :, :-1] + by[:, :, 1:])

    return (q0, bx, by)
import numpy as np
import scipy.sparse as sp
from dataclasses import replace
from .initialize import initialize_1d, initialize_2d, VariableState, ModelState

# =========================
# Downsample operators
# =========================

def downsample_phi_2d(v):
    Mx = v.shape[0] - 1
    My = v.shape[1] - 1

    Mxc = Mx // 2
    Myc = My // 2

    vc = np.zeros((Mxc + 1, Myc + 1), order='F')

    ind_center = slice(2, Mx, 2)
    ind_prev = slice(1, Mx - 1, 2)
    ind_next = slice(3, Mx + 1, 2)

    v_cc = v[ind_center, ind_center]  # center-center

    v_cl = v[ind_center, ind_prev]  # center-left
    v_cr = v[ind_center, ind_next]  # center-right
    v_uc = v[ind_prev, ind_center]  # up-center
    v_dc = v[ind_next, ind_center]  # down-center

    v_ul = v[ind_prev, ind_prev]  # up-left
    v_ur = v[ind_prev, ind_next]  # up-right
    v_dl = v[ind_next, ind_prev]  # down-left
    v_dr = v[ind_next, ind_next]  # down-right

    vc[1:Mxc, 1:Myc] = (4 * v_cc +
                        2 * (v_cl + v_cr + v_uc + v_dc) +
                        (v_ul + v_ur + v_dl + v_dr)) / 16.0

    vc[0, 1:Myc] = (4 * v[0, ind_center] +
                    2 * (v[1, ind_center] + v[0, ind_prev] + v[0, ind_next]) +
                    (v[1, ind_prev] + v[1, ind_next])) / 12.0

    vc[Mxc, 1:Myc] = (4 * v[Mx, ind_center] +
                      2 * (v[Mx - 1, ind_center] + v[Mx, ind_prev] + v[Mx, ind_next]) +
                      (v[Mx - 1, ind_prev] + v[Mx - 1, ind_next])) / 12.0

    vc[1:Mxc, 0] = (4 * v[ind_center, 0] +
                    2 * (v[ind_prev, 0] + v[ind_next, 0] + v[ind_center, 1]) +
                    (v[ind_prev, 1] + v[ind_next, 1])) / 12.0

    vc[1:Mxc, Myc] = (4 * v[ind_center, My] +
                      2 * (v[ind_prev, My] + v[ind_next, My] + v[ind_center, My - 1]) +
                      (v[ind_prev, My - 1] + v[ind_next, My - 1])) / 12.0

    vc[0, 0] = (4 * v[0, 0] + 2 * (v[1, 0] + v[0, 1]) + v[1, 1]) / 9.0
    vc[0, Myc] = (4 * v[0, My] + 2 * (v[1, My] + v[0, My - 1]) + v[1, My - 1]) / 9.0
    vc[Mxc, 0] = (4 * v[Mx, 0] + 2 * (v[Mx - 1, 0] + v[Mx, 1]) + v[Mx - 1, 1]) / 9.0
    vc[Mxc, Myc] = (4 * v[Mx, My] + 2 * (v[Mx - 1, My] + v[Mx, My - 1]) + v[Mx - 1, My - 1]) / 9.0

    return vc


def downsample_phi_1d(v):
    phi = np.asarray(v, dtype=np.float64).reshape(-1, order='F')
    if phi.size < 2:
        raise ValueError("phi must contain at least two samples")

    M = phi.shape[0] - 1
    if M % 2 != 0:
        raise ValueError("phi length minus one must be even for downsampling")

    Mc = M // 2

    phic = np.zeros((Mc + 1, 1), dtype=np.float64, order='F')

    if Mc > 0:
        phic[0, 0] = (2.0 / 3.0) * phi[0] + (1.0 / 3.0) * phi[1]
        phic[-1, 0] = (1.0 / 3.0) * phi[-2] + (2.0 / 3.0) * phi[-1]

    if Mc > 1:
        phic[1:Mc, 0] = (
            0.5 * phi[2:-1:2]
            + 0.25 * (phi[1:-2:2] + phi[3::2])
        )

    return phic


def _prolong_mat_linear(nc: int) -> sp.csr_matrix:
    if nc < 1:
        raise ValueError("linear prolongation requires at least one coarse point")

    nr = 2 * (nc - 1) + 1
    rows = [np.arange(0, nr, 2)]
    cols = [np.arange(nc)]
    data = [np.ones(nc, dtype=np.float64)]

    if nc > 1:
        rows_even = np.arange(1, nr - 1, 2)
        rows.append(np.repeat(rows_even, 2))
        cols_even = np.empty(rows_even.size * 2, dtype=np.int64)
        cols_even[0::2] = np.arange(nc - 1)
        cols_even[1::2] = np.arange(1, nc)
        cols.append(cols_even)
        data.append(np.full(rows_even.size * 2, 0.5, dtype=np.float64))

    rows_all = np.concatenate(rows)
    cols_all = np.concatenate(cols)
    data_all = np.concatenate(data)

    return sp.coo_matrix((data_all, (rows_all, cols_all)), shape=(nr, nc)).tocsr()


def _prolong_mat_nearest(nc: int) -> sp.csr_matrix:
    if nc < 0:
        raise ValueError("nearest prolongation size must be non-negative")

    nr = 2 * nc
    if nr == 0:
        return sp.csr_matrix((0, nc), dtype=np.float64)

    rows = np.arange(nr, dtype=np.int64)
    cols = np.repeat(np.arange(nc, dtype=np.int64), 2)
    data = np.ones(nr, dtype=np.float64)

    return sp.coo_matrix((data, (rows, cols)), shape=(nr, nc)).tocsr()


def _restriction_from_prolong(prolong: sp.csr_matrix) -> sp.csr_matrix:
    col_sums = np.asarray(prolong.sum(axis=0)).ravel()
    if col_sums.size and np.any(col_sums == 0):
        raise ValueError("prolongation matrix contains zero-sum columns")

    if col_sums.size:
        scale = sp.diags(1.0 / col_sums)
    else:
        scale = sp.csr_matrix((prolong.shape[1], prolong.shape[1]), dtype=np.float64)
    return (prolong @ scale).transpose().tocsr()


def downsample_q_2d(q: np.ndarray, nt: int, nx: int, ny: int) -> np.ndarray:
    q_vec = np.asarray(q, dtype=np.float64).reshape(-1, 1)

    len_t = (nt - 1) * nx * ny
    len_x = nt * (nx - 1) * ny
    len_y = nt * nx * (ny - 1)
    expected = len_t + len_x + len_y
    if q_vec.shape[0] != expected:
        raise ValueError(f"q length {q_vec.shape[0]} does not match expected size {expected}")

    nt2 = (nt + 1) // 2
    nx2 = (nx + 1) // 2
    ny2 = (ny + 1) // 2

    prolong_t = sp.kron(
        sp.kron(_prolong_mat_nearest(nt2 - 1), _prolong_mat_linear(nx2), format='csr'),
        _prolong_mat_linear(ny2),
        format='csr'
    )
    prolong_x = sp.kron(
        sp.kron(_prolong_mat_linear(nt2), _prolong_mat_nearest(nx2 - 1), format='csr'),
        _prolong_mat_linear(ny2),
        format='csr'
    )
    prolong_y = sp.kron(
        sp.kron(_prolong_mat_linear(nt2), _prolong_mat_linear(nx2), format='csr'),
        _prolong_mat_nearest(ny2 - 1),
        format='csr'
    )

    restri_t = _restriction_from_prolong(prolong_t)
    restri_x = _restriction_from_prolong(prolong_x)
    restri_y = _restriction_from_prolong(prolong_y)

    qt = restri_t @ q_vec[:len_t]
    qx = restri_x @ q_vec[len_t:len_t + len_x]
    qy = restri_y @ q_vec[len_t + len_x:]

    return np.vstack((qt, qx, qy))


# =========================
# Interpolate operators of 1D
# =========================

def _interpolate_linear_1d(f, axis):
    slices_pre = [slice(None)] * f.ndim
    slices_next = [slice(None)] * f.ndim
    slices_pre[axis] = slice(0, -1)
    slices_next[axis] = slice(1, None)

    return 0.5 * (f[tuple(slices_pre)] + f[tuple(slices_next)])


def _interpolate_t_stagger_1d(f):
    nx, nt = f.shape

    nxR = 2 * (nx - 1) + 1
    ntR = 2 * nt

    fR = np.zeros((nxR, ntR), dtype=np.float64, order='F')
    fR[0::2, 0:-1:2] = f
    fR[0::2, 1::2] = f

    if nxR > 1:
        fR[1:-1:2, :] = _interpolate_linear_1d(fR[0::2, :], 0)

    return fR


def _interpolate_phi_1d(phi_vec, nx, nt):
    phi = phi_vec.reshape((nx, nt), order='F')

    nxR = 2 * (nx - 1) + 1
    ntR = 2 * (nt - 1) + 1

    phiR = np.zeros((nxR, ntR), dtype=np.float64, order='F')

    phiR[0::2, 0::2] = phi
    if nxR > 1:
        phiR[1:-1:2, 0::2] = _interpolate_linear_1d(phiR[0::2, 0::2], 0)
    if ntR > 1:
        phiR[:, 1:-1:2] = _interpolate_linear_1d(phiR[:, 0::2], 1)

    return phiR.flatten(order='F').reshape(-1, 1)


def _interpolate_z_1d(z_mat, nx, nt):
    nt_z = nt - 1

    nxR = 2 * (nx - 1) + 1
    ntR_z = 2 * nt_z

    len_fine = nxR * ntR_z
    num_cols = z_mat.shape[1]
    zR = np.zeros((len_fine, num_cols), dtype=np.float64, order='F')

    for j in range(num_cols):
        col = z_mat[:, j].reshape((nx, nt_z), order='F')
        interp = _interpolate_t_stagger_1d(col)
        zR[:, j] = interp.flatten(order='F')

    return zR


def _interpolate_var_1d(var: VariableState, model: ModelState):
    nx = model.nx
    nt = model.nt

    phi_interp = _interpolate_phi_1d(var.phi, nx, nt)
    beta_interp = _interpolate_z_1d(var.beta, nx, nt)

    return replace(var, phi=phi_interp, beta=beta_interp)





# =========================
# Interpolate operators of 2D
# =========================
def _interpolate_linear_2d(f, axis):
    slices_pre = [slice(None)] * f.ndim
    slices_next = [slice(None)] * f.ndim
    slices_pre[axis] = slice(0, -1)
    slices_next[axis] = slice(1, None)

    return 0.5 * (f[tuple(slices_pre)] + f[tuple(slices_next)])


def _interpolate_t_stagger_2d(f):
    ny, nx, nt = f.shape

    nyR = 2 * (ny - 1) + 1
    nxR = 2 * (nx - 1) + 1
    ntR = 2 * nt

    fR = np.zeros((nyR, nxR, ntR), dtype=np.float64, order='F')
    fR[0::2, 0::2, 0:-1:2] = f  # t - nearest
    fR[0::2, 0::2, 1::2] = f  # t - nearest (evenT)
    fR[1:-1:2, 0::2, :] = _interpolate_linear_2d(fR[0::2, 0::2, :], 0)
    fR[:, 1:-1:2, :] = _interpolate_linear_2d(fR[:, 0::2, :], 1)

    return fR


def _interpolate_phi_2d(phi_vec, ny, nx, nt):
    phi = phi_vec.reshape((ny, nx, nt), order='F')

    nyR = 2 * (ny - 1) + 1
    nxR = 2 * (nx - 1) + 1
    ntR = 2 * (nt - 1) + 1

    phiR = np.zeros((nyR, nxR, ntR), dtype=np.float64, order='F')

    phiR[0::2, 0::2, 0::2] = phi
    phiR[1:-1:2, 0::2, 0::2] = _interpolate_linear_2d(phiR[0::2, 0::2, 0::2], 0)
    phiR[:, 1:-1:2, 0::2] = _interpolate_linear_2d(phiR[:, 0::2, 0::2], 1)
    phiR[:, :, 1:-1:2] = _interpolate_linear_2d(phiR[:, :, 0::2], 2)

    return phiR.flatten(order='F').reshape(-1, 1)


def _interpolate_z_2d(z_mat, ny, nx, nt):
    nt_z = nt - 1

    nyR = 2 * (ny - 1) + 1
    nxR = 2 * (nx - 1) + 1
    ntR_z = 2 * nt_z

    len_fine = nyR * nxR * ntR_z
    zR = np.zeros((len_fine, 10), dtype=np.float64, order='F')

    for j in range(10):
        col = z_mat[:, j]
        col_3d = col.reshape((ny, nx, nt_z), order='F')
        interp_3d = _interpolate_t_stagger_2d(col_3d)
        zR[:, j] = interp_3d.flatten(order='F')

    return zR


def _interpolate_var_2d(var: VariableState, model: ModelState):
    nx = model.nx
    ny = model.ny
    nt = model.nt

    phi_interp = _interpolate_phi_2d(var.phi, ny, nx, nt)
    beta_interp = _interpolate_z_2d(var.beta, ny, nx, nt)

    return replace(var, phi=phi_interp, beta=beta_interp)


# =========================
# Interface - Jump to the next level
# =========================

def jump_next_level(var: VariableState, model: ModelState, rho0, rho1, nt_fine, weight=None) -> tuple[VariableState, ModelState]:
    dim = model.dim

    if dim == 1:
        varR = _interpolate_var_1d(var, model)
        var_init, modelR = initialize_1d(rho0, rho1, nt_fine, weight=weight)
    elif dim == 2:
        varR = _interpolate_var_2d(var, model)
        var_init, modelR = initialize_2d(rho0, rho1, nt_fine, weight=weight)
    else:
        raise NotImplementedError(f"Dimension {dim} not supported for multilevel interpolation.")

    varR.qInd = var_init.qInd.copy()
    varR.z = var_init.z.copy()
    varR.q = modelR.grad @ varR.phi
    varR.alpha = var_init.alpha.copy()

    if weight is not None:
        varR.q /= weight.reshape(-1, 1)

    return varR, modelR

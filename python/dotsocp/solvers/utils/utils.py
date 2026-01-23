import numpy as np


def integral_l2(f, h=None):
    f = f.reshape(-1, 1, order='F')
    if h is None:
        n = f.shape[0]
        h = 1.0 / n
    return h * np.sum(f)


def check_mass_conservation(rho, tol=1e-2, show=False):
    ny, nx, nt = rho.shape

    rho2 = rho.reshape(nx * ny, nt, order='F')
    sum_rho = np.sum(rho2, axis=0)

    nega_rho = rho2.copy()
    nega_rho[rho2 >= 0] = 0
    sum_nega_rho = np.sum(nega_rho, axis=0)

    err = np.max(np.abs(np.concatenate([sum_rho - 1, sum_nega_rho])))

    flag = 1 if err <= tol else 0

    if show:
        print("Total mass (sum of rho):")
        print(sum_rho.flatten()[:10])
        if flag == 0:
            print(f"Warning: Mass conservation violation exceeds {tol}")

    return flag

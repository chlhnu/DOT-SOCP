from dataclasses import dataclass
from typing import Dict, Optional, Protocol
from numpy import ndarray

@dataclass
class SolverParams:
    tau: float
    sigma: float
    maxit: int
    tol: float
    tol_feasorg: float
    time_limit: float
    use_scaling: bool


@dataclass
class ScalingControl:
    level: int
    first_iter: int
    second_iter: int
    check_iters: int
    ratio_threshold: float
    max_feas: float
    rel_gap: float
    sigma_scale: float


class CppExtProtocol(Protocol):
    """Wrapper for C++ extension module with bound methods."""
    def oper_bfd(
        self,
        z, q,
        scale_bf: float = ...,
        scale_d: float = ...
    ) -> None: ...
    # In-place modification of z given q

    def oper_bfd_conj(
        self,
        q, z,
        scale_bf: float = ...
    ) -> None: ...
    # In-place modification of q given z

    def proj_soc(self, z, zz) -> None: ... 
    # In-place modification of z given zz


@dataclass
class ProblemContext:
    dim: int
    h: float
    qInd: Dict[str, int]
    At: ndarray
    grad: ndarray
    c: ndarray
    kernel: ndarray
    diag_q_inv: ndarray
    scale_bf: float
    D: float
    E: float
    norm_c_ref: float
    norm_d_ref: float
    cpp_ext: CppExtProtocol  # Placeholder for C++ extension module
    nt: int
    nx: int
    ny: Optional[int] = None  # invalid for 1D problems
    weight: Optional[ndarray] = None  # for uneven OT problems


@dataclass
class IterationState:
    phi: ndarray
    q: ndarray
    z: ndarray
    alpha: ndarray
    beta: ndarray
    z2: ndarray
    q2: ndarray
    proj_z_beta: ndarray
    c_scaled: ndarray
    c_scale: float
    d_scale: float
    tmp_q: ndarray
    resi_alpha: ndarray
    resi_beta: ndarray


@dataclass
class KKtResult:
    org_kkt_resi: list
    org_pd_ratio: float
    org_max_kkt: float
    kkt_resi: list
    pd_ratio: float
    max_kkt: float
    pd_gap: float
    validation: float

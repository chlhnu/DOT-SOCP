import numpy as np
import os

from .utils.metrics import norm_l2, fnorm_l2
from .utils.dataclass import (
    ProblemContext,
    IterationState,
    KKtResult,
)
from .utils.admm_utils import RunHistoryManager


def _compute_org_dot_terms_2d(
    state: IterationState,
    context: ProblemContext,
    dual_scale: float,
):
    """Compute rhoT, rhoFq, rhoBx, rhoBy terms for 2D DOT KKT metrics."""
    qInd = context.qInd
    D, E = context.D, context.E
    weight = 1.0 if context.weight is None else context.weight
    ny = context.ny
    nx = context.nx
    nt = context.nt

    assert ny is not None, "ny must be specified for 2D problems."

    weight_head = 1.0 if isinstance(weight, float) else weight[:qInd['bx']]
    alpha_head = state.alpha[:qInd['bx']]
    q_head = state.q[:qInd['bx']]
    z2_bulk = state.z2[:, 1:-1]

    weight_arr = np.asarray(weight_head, dtype=np.float64).reshape(-1, 1, order="F")
    alpha_arr = np.asarray(alpha_head, dtype=np.float64).reshape(-1, 1, order="F")
    q_arr     = np.asarray(q_head,   dtype=np.float64).reshape(-1, 1, order="F")
    z_arr     = np.asarray(z2_bulk,  dtype=np.float64)

    rho_t    = (dual_scale * state.c_scale * D) * (weight_arr * alpha_arr)
    scaled_q = (state.d_scale / D) * q_arr
    scaled_z = (state.d_scale / E) * z_arr

    quadratic = np.sum(scaled_z ** 2, axis=1, keepdims=True) / 4.0
    rho_fq = np.maximum(rho_t + scaled_q + quadratic, 0.0)

    rhoT_3d = rho_t.reshape((ny, nx, nt - 1), order='F')
    zeros_pad = np.zeros((ny, nx, 1), order='F')
    rho_padded = np.concatenate((zeros_pad, rhoT_3d, zeros_pad), axis=2)
    rho = (rho_padded[:, :, :-1] + rho_padded[:, :, 1:]) / 2.0

    rho_avg_x = (rho[:, :-1, :] + rho[:, 1:, :]) / 2.0
    rho_avg_x_flat = rho_avg_x.reshape(-1, 1, order='F')
    q_bx = state.q[qInd['bx']: qInd['by']]
    rhoBx = (state.d_scale / D) * (rho_avg_x_flat * q_bx)

    rho_avg_y = (rho[:-1, :, :] + rho[1:, :, :]) / 2.0
    rho_avg_y_flat = rho_avg_y.reshape(-1, 1, order='F')
    q_by = state.q[qInd['by']:]
    rhoBy = (state.d_scale / D) * (rho_avg_y_flat * q_by)

    return rho_t, rho_fq, rhoBx, rhoBy


def _compute_org_dot_terms_1d(
    state: IterationState,
    context: ProblemContext,
    dual_scale: float,
):
    """Compute rhoT, rhoFq, rhoBx, and mx terms for 1D DOT KKT metrics."""
    qInd = context.qInd
    D, E = context.D, context.E
    weight = 1.0 if context.weight is None else context.weight
    nx = context.nx
    nt = context.nt

    weight_head = 1.0 if isinstance(weight, float) else weight[:qInd['bx']]
    alpha_head = state.alpha[:qInd['bx']]
    q_head = state.q[:qInd['bx']]
    z2_bulk = state.z2[:, 1:5]

    weight_arr = np.asarray(weight_head, dtype=np.float64).reshape(-1, 1, order="F")
    alpha_arr = np.asarray(alpha_head, dtype=np.float64).reshape(-1, 1, order='F')
    q_arr = np.asarray(q_head, dtype=np.float64).reshape(-1, 1, order='F')
    z_arr = np.asarray(z2_bulk, dtype=np.float64)

    rho_t    = (dual_scale * state.c_scale * D) * (weight_arr * alpha_arr)
    scaled_q = (state.d_scale / D) * q_arr
    scaled_z = (state.d_scale / E) * z_arr
    quadratic = np.sum(scaled_z ** 2, axis=1, keepdims=True) / 4.0
    rho_fq = np.maximum(rho_t + scaled_q + quadratic, 0.0)

    rhoT_2d = rho_t.reshape((nx, nt - 1), order='F')
    zeros_pad = np.zeros((nx, 1), order='F')
    rho_padded = np.concatenate((zeros_pad, rhoT_2d, zeros_pad), axis=1)
    rho = (rho_padded[:, :-1] + rho_padded[:, 1:]) / 2.0

    rho_avg_x = (rho[:-1, :] + rho[1:, :]) / 2.0
    rho_avg_x_flat = rho_avg_x.reshape(-1, 1, order='F')
    q_tail = state.q[qInd['bx']:]
    rhoBx = (state.d_scale / D) * (rho_avg_x_flat * q_tail)

    return rho_t, rho_fq, rhoBx


def kkt_metric_dot(
    state: IterationState,
    context: ProblemContext,
    dual_scale: float = 1.0,
) -> KKtResult:
    qInd = context.qInd
    is_check_rho_b = os.getenv("DOTSOCP_CHECK_RHO_B", "1") == "1"

    context.cpp_ext.oper_bfd_conj(state.q2, state.beta, context.scale_bf)
    context.cpp_ext.proj_soc(state.proj_z_beta, state.z - dual_scale * state.beta)
    weighted_alpha = state.alpha if context.weight is None else context.weight * state.alpha

    norm_q_val = norm_l2(state.q, context.h)
    norm_z_val = fnorm_l2(state.z, context.h)
    norm_aphi = norm_l2(state.tmp_q, context.h)
    norm_alpha_real = dual_scale * norm_l2(state.alpha, context.h)
    norm_beta_real = dual_scale * fnorm_l2(state.beta, context.h)
    norm_FBbeta = dual_scale * norm_l2(state.q2, context.h)

    prim_fea1 = norm_l2(state.resi_alpha, context.h)
    prim_fea2 = fnorm_l2(state.resi_beta, context.h)
    dual_fea1 = dual_scale * norm_l2(context.At @ state.alpha - state.c_scaled, context.h)
    complem = fnorm_l2(state.z - state.proj_z_beta, context.h)
    dual_fea2 = dual_scale * norm_l2(state.q2 + weighted_alpha, context.h)

    if context.dim == 1:
        rhoT, rhoFq, rhoBx = _compute_org_dot_terms_1d(state, context, dual_scale)
        normRho = norm_l2(rhoT, context.h)
        norm_rhoFq = norm_l2(rhoFq, context.h)

        mx = (dual_scale * state.c_scale * context.D) * weighted_alpha[qInd['bx']:]
        normM = norm_l2(mx, context.h)
        norm_rho_b = norm_l2(rhoBx, context.h)

        dotcomple = norm_l2(rhoT - rhoFq, context.h)
        mRhoB = norm_l2(mx - rhoBx, context.h)
    elif context.dim == 2:
        rhoT, rhoFq, rhoBx, rhoBy = _compute_org_dot_terms_2d(state, context, dual_scale)
        normRho = norm_l2(rhoT, context.h)
        norm_rhoFq = norm_l2(rhoFq, context.h)

        mx = (dual_scale * state.c_scale * context.D) * weighted_alpha[qInd['bx']: qInd['by']]
        my = (dual_scale * state.c_scale * context.D) * weighted_alpha[qInd['by']:]
        normM = np.sqrt(norm_l2(mx, context.h) ** 2 + norm_l2(my, context.h) ** 2)
        norm_rho_b = np.sqrt(norm_l2(rhoBx, context.h) ** 2 + norm_l2(rhoBy, context.h) ** 2)

        dotcomple = norm_l2(rhoT - rhoFq, context.h)
        mRhoB = np.sqrt(norm_l2(mx - rhoBx, context.h) ** 2 + norm_l2(my - rhoBy, context.h) ** 2)
    else:
        raise NotImplementedError(f"Dimension {context.dim} not supported for KKT metrics.")

    kkt_const = 1.0

    # KKT residuals of original problem
    kkt_resi_org = [
        prim_fea1 / (kkt_const * context.D / state.d_scale + norm_aphi + norm_q_val),
        prim_fea2 / (kkt_const * context.E / state.d_scale + context.norm_d_ref),
        dual_fea1 / (kkt_const / state.c_scale + context.norm_c_ref),
        complem / (kkt_const * context.E / state.d_scale + norm_z_val + norm_beta_real),
        dual_fea2 / (kkt_const + norm_FBbeta + norm_alpha_real),
        dotcomple / (kkt_const + normRho + norm_rhoFq),
        mRhoB / (kkt_const + normM + norm_rho_b),
    ]
    pd_ratio_org = max(kkt_resi_org[0], kkt_resi_org[1]) / (max(kkt_resi_org[2], kkt_resi_org[4]) + 1e-16)

    # KKT residuals of scaled problem
    kkt_resi = [
        prim_fea1 / (kkt_const + norm_aphi + norm_q_val),
        prim_fea2 / (kkt_const + context.norm_d_ref),
        dual_fea1 / (kkt_const + context.norm_c_ref),
        complem / (kkt_const + norm_z_val + norm_beta_real),
        dual_fea2 / (kkt_const + norm_FBbeta + norm_alpha_real),
    ]
    pd_ratio = max(kkt_resi[0], kkt_resi[1]) / (max(kkt_resi[2], kkt_resi[4]) + 1e-16)

    # Objective
    prefactor = dual_scale * state.c_scale * state.d_scale * context.h
    weighted_q = state.q if context.weight is None else context.weight * state.q
    pri_val = prefactor * np.vdot(weighted_q, state.alpha)
    dual_val = prefactor * np.vdot(state.c_scaled, state.phi)
    pd_gap = np.abs(pri_val - dual_val) / (1 + np.abs(pri_val) + np.abs(dual_val))

    # The kkt errors used to validate convergence
    if is_check_rho_b:
        validated_indices = [0, 2, 5, 6]
    else:
        validated_indices = [0, 2, 5]
    
    validation = max(kkt_resi_org[i] for i in validated_indices)

    return KKtResult(
        org_kkt_resi=kkt_resi_org,
        org_pd_ratio=pd_ratio_org,
        org_max_kkt=max(kkt_resi_org),
        kkt_resi=kkt_resi,
        pd_ratio=pd_ratio,
        max_kkt=max(kkt_resi),
        pd_gap=pd_gap,
        validation=validation,
    )


def make_run_history_manager(max_records: int) -> RunHistoryManager:
    is_check_rho_b = os.getenv("DOTSOCP_CHECK_RHO_B", "1") == "1"

    # The kkt errors used to validate convergence
    if is_check_rho_b:
        validated_indices = [0, 2, 5, 6]
    else:
        validated_indices = [0, 2, 5]
    
    kkt_labels = [
        "Prim feas alpha",
        "Prim feas beta",
        "Dual feas c",
        "SOCP compl",
        "Dual feas q",
        "DOT compl",
        "Momentum rhoB",
    ]
    kkt_short = ["Prim-a", "Prim-b", "Dual-c", "SOCP-compl", "Dual-q", "DOT-compl", "rhoB"]

    kkt_labels = [label + " (*)" if idx in validated_indices else label for idx, label in enumerate(kkt_labels)]
    kkt_short = [label + "(*)" if idx in validated_indices else label for idx, label in enumerate(kkt_short)]

    manager = RunHistoryManager(
        max_record_numbers=max_records,
        kkt_labels=kkt_labels,
        name="KKT errors for DOT-SOCP",
        kkt_short_labels=kkt_short,
    )
    return manager

import time
import logging
import copy
from typing import Dict, Callable
import numpy as np

from .algo_updates import inexact_palm_update
from .kkt_metrics import kkt_metric_dot, make_run_history_manager
from .utils.metrics import norm_l2, fnorm_l2
from .utils.admm_utils import PenaltyParamManager
from .utils.backend import make_socp_backend
from .utils.operators import (
    oper_q_1d, create_fft_kernel_2d,
    oper_q_2d, create_fft_kernel_3d,
)
from .utils.dataclass import (
    SolverParams,
    ScalingControl,
    ProblemContext,
    IterationState,
)
from ..utils.initialize import VariableState, ModelState


def _build_solver_params(opts: Dict) -> SolverParams:
    tau = opts.get('tau', 1.9)
    sigma = opts.get('sigma', 1.0)
    maxit = opts.get('maxit', 3000)
    tol = opts.get('tol', 1e-4)
    tol_feasorg = 5 * tol
    time_limit = opts.get('time_limit', 3600)
    use_scaling = opts.get('scaling', True)
    return SolverParams(
        tau=tau,
        sigma=sigma,
        maxit=maxit,
        tol=tol,
        tol_feasorg=tol_feasorg,
        time_limit=time_limit if time_limit > 0 else float('inf'),
        use_scaling=use_scaling,
    )


def _init_scaling_control(use_scaling: bool) -> ScalingControl:
    if use_scaling:
        return ScalingControl(
            level=1,
            first_iter=10,
            second_iter=50,
            check_iters=100,
            ratio_threshold=1.2,
            max_feas=float('inf'),
            rel_gap=float('inf'),
            sigma_scale=1.0,
        )
    return ScalingControl(
        level=0,
        first_iter=10,
        second_iter=50,
        check_iters=100,
        ratio_threshold=1.2,
        max_feas=float('inf'),
        rel_gap=float('inf'),
        sigma_scale=1.0,
    )


def _build_problem_context(var: VariableState, model: ModelState) -> ProblemContext:
    D = getattr(var, 'D', 1.0)
    E = getattr(var, 'E', 1.0)
    dim = model.dim
    weight = getattr(model, 'weight', None)

    if dim == 2:
        nt, nx, ny = model.nt, model.nx, model.ny
        assert ny is not None, "ny must be specified for 2D problems."

        h = 1.0 / (nx * ny * nt)
        kernel = (D ** 2) * create_fft_kernel_3d(nt, nx, ny)
        diag_q = oper_q_2d(ny, nx, nt, D, E, weight=weight)
        diag_q_inv = (1.0 / diag_q).reshape(-1, 1, order='F')

        backend = make_socp_backend(dim=2, nt=nt, nx=nx, ny=ny)

    elif dim == 1:
        nt, nx = model.nt, model.nx
        ny = None
        h = 1.0 / (nx * nt)
        kernel = (D ** 2) * create_fft_kernel_2d(nt, nx)
        diag_q = oper_q_1d(nx, nt, D, E, weight=weight)
        diag_q_inv = (1.0 / diag_q).reshape(-1, 1, order='F')

        backend = make_socp_backend(dim=1, nt=nt, nx=nx)
    
    norm_c_ref = norm_l2(model.c, h) * np.sqrt(nt)
    norm_d_ref = 1.0

    return ProblemContext(
        dim=model.dim,
        backend=backend,
        nx=nx,
        ny=ny,
        nt=nt,
        h=h,
        qInd=var.qInd,
        At=model.grad.T,
        grad=model.grad,
        c=model.c,
        kernel=kernel,
        diag_q_inv=diag_q_inv,
        scale_bf=E / D,
        D=D,
        E=E,
        norm_c_ref=norm_c_ref,
        norm_d_ref=norm_d_ref,
        weight=weight,
    )


def _initialize_iteration_state(var: VariableState, sigma, context: ProblemContext) -> IterationState:
    # Initialize variables
    phi = var.phi.copy(order='F')
    q = var.q.copy(order='F')
    z = var.z.copy(order='F')
    alpha = var.alpha.copy(order='F') / sigma
    beta = var.beta.copy(order='F') / sigma
    c_scale = getattr(var, 'cScale', 1.0)
    d_scale = getattr(var, 'dScale', 1.0)
    c_scaled = np.asfortranarray(context.c / sigma)

    if context.dim == 1:
        len_z = (context.nt - 1) * context.nx
        z2 = np.zeros((len_z, 6), dtype=np.float64, order='F')
    elif context.dim == 2:
        assert context.ny is not None, "ny must be specified for 2D problems."
        len_z = (context.nt - 1) * context.nx * context.ny
        z2 = np.zeros((len_z, 10), dtype=np.float64, order='F')
    else:
        raise NotImplementedError(f"Dimension {context.dim} not supported for iteration state initialization.")

    q2 = np.zeros_like(q, order='F')
    proj_z_beta = np.zeros_like(z, order='F')
    # tmp_q = np.zeros_like(q, order='F')
    resi_alpha = np.zeros_like(q, order='F')
    resi_beta = np.zeros_like(z, order='F')

    # Precompute
    tmp_q = context.grad @ phi
    context.backend.oper_bfd(z, tmp_q, context.scale_bf, context.E / d_scale)
    context.backend.oper_bfd(z2, q, context.scale_bf, context.E / d_scale)
    context.backend.oper_bfd_conj(alpha, -beta, context.scale_bf)

    if context.weight is not None:
        alpha /= context.weight

    return IterationState(
        phi=phi,
        q=q,
        z=z,
        alpha=alpha,
        beta=beta,
        z2=z2,
        q2=q2,
        proj_z_beta=proj_z_beta,
        c_scaled=c_scaled,
        c_scale=c_scale,
        d_scale=d_scale,
        tmp_q=tmp_q,
        resi_alpha=resi_alpha,
        resi_beta=resi_beta,
    )


def _compute_rescale_norms(state: IterationState, context: ProblemContext, sigma: float) -> Dict[str, float]:
    norm_phi = norm_l2(state.phi, context.h)
    norm_q_val = norm_l2(state.q, context.h)
    norm_z_val = fnorm_l2(state.z, context.h)
    norm_alpha_val = sigma * norm_l2(state.alpha, context.h)
    norm_beta_val = sigma * fnorm_l2(state.beta, context.h)
    norm_phis = max(norm_phi, norm_q_val, norm_z_val)
    norm_alps = max(norm_alpha_val, norm_beta_val)
    min_norm = min(norm_alps, norm_phis)
    ratio = max(norm_alps, norm_phis) / (min_norm + 1e-16)
    return {
        'norm_phi': norm_phi,
        'norm_q': norm_q_val,
        'norm_z': norm_z_val,
        'norm_alpha': norm_alpha_val,
        'norm_beta': norm_beta_val,
        'norm_phis': norm_phis,
        'norm_alps': norm_alps,
        'ratio': ratio,
    }


def _apply_rescale(
    state: IterationState,
    context: ProblemContext,
    sigma: float,
    norms: Dict[str, float],
    scaling_ctrl: ScalingControl,
    iter_idx: int,
):
    d_scale2 = norms['norm_phis']
    c_scale2 = norms['norm_alps']

    scale_ratio = c_scale2 / (d_scale2 + 1e-16)
    c_factor = d_scale2 / ((c_scale2 ** 2) + 1e-16)

    state.c_scaled = np.asfortranarray(state.c_scaled * c_factor)
    context.norm_c_ref /= c_scale2
    context.norm_d_ref /= d_scale2

    state.phi = np.asfortranarray(state.phi / d_scale2)
    state.tmp_q = np.asfortranarray(state.tmp_q / d_scale2) # tmp_q == grad @ phi
    state.q = np.asfortranarray(state.q / d_scale2)
    state.z = np.asfortranarray(state.z / d_scale2)
    state.alpha = np.asfortranarray(state.alpha * c_factor)
    state.beta = np.asfortranarray(state.beta * c_factor)

    sigma *= scale_ratio
    scaling_ctrl.sigma_scale *= scale_ratio

    state.d_scale *= d_scale2
    state.c_scale *= c_scale2
    scale_d = context.E / state.d_scale
    context.backend.oper_bfd(state.z2, state.q, context.scale_bf, scale_d)

    logging.debug(f"   [Rescale #{scaling_ctrl.level - 1}] Iter: {iter_idx}, Sigma: {sigma:.2e}, Ratio: {scale_ratio:.2f}")
    scaling_ctrl.level += 1

    return sigma


def _maybe_rescale(
    iter_idx: int,
    sigma: float,
    state: IterationState,
    context: ProblemContext,
    scaling_ctrl: ScalingControl,
):
    if scaling_ctrl.level == 0:
        return sigma
    norms = None
    scale_yes = False
    if scaling_ctrl.level >= 3 and (iter_idx % scaling_ctrl.check_iters == 0):
        norms = _compute_rescale_norms(state, context, sigma)
        if norms['ratio'] > scaling_ctrl.ratio_threshold:
            scale_yes = True
    condition_1 = (
        scaling_ctrl.level == 1
        and scaling_ctrl.max_feas < 2e-2
        and iter_idx >= scaling_ctrl.first_iter
        and scaling_ctrl.rel_gap < 5e-2
    )
    condition_2 = (
        scaling_ctrl.level == 2
        and scaling_ctrl.max_feas < 5e-3
        and iter_idx >= scaling_ctrl.second_iter
        and scaling_ctrl.rel_gap < 1e-2
    )
    if (condition_1 or condition_2) and not scale_yes:
        norms = _compute_rescale_norms(state, context, sigma)
        scale_yes = True
    if scale_yes:
        sigma = _apply_rescale(
            state,
            context,
            sigma,
            norms or _compute_rescale_norms(state, context, sigma),
            scaling_ctrl,
            iter_idx,
        )
    return sigma


def run_socp_solver(
        var: VariableState, 
        model: ModelState,
        opts,
        update_fn: Callable = inexact_palm_update,
        kkt_fn: Callable = kkt_metric_dot,
    ):
    params = _build_solver_params(opts)
    scaling_ctrl = _init_scaling_control(params.use_scaling)
    context = _build_problem_context(var, model)

    penalty_manager = PenaltyParamManager()
    history_manager = make_run_history_manager(params.maxit)
    history_manager.start()
    history_manager.create_tol_progress(target_tol=params.tol)

    use_feas_org = False
    is_converged = False
    elapsed = 0.0

    total_start_time = time.perf_counter()
    sigma = params.sigma
    state = _initialize_iteration_state(var, sigma, context)

    # Main iteration loop
    for it in range(1, params.maxit + 1):
        sigma = _maybe_rescale(
            it,
            sigma,
            state,
            context,
            scaling_ctrl,
        )

        update_fn(state, context, params, history_manager)

        adjust_sigma_yes = penalty_manager.is_to_adjust(it)
        needs_kkt = adjust_sigma_yes or (it == params.maxit)

        if needs_kkt:
            kkt_result = kkt_fn(
                state=state,
                context=context,
                dual_scale=sigma,
            )
            history_manager.record(
                current_it=it,
                kkt_errors=kkt_result.org_kkt_resi,
                history={
                    'pd_gap': kkt_result.pd_gap,
                },
            )
            history_manager.show_tol_progress(it, kkt_result.validation)

            # logging.info(
            #     f"Iter {it}, KKT Error {kkt_result.max_kkt:.2e}, Primal/Dual Ratio {kkt_result.pd_ratio:.2f}"
            # )

            if kkt_result.validation < params.tol:
                is_converged = True
                break

            if scaling_ctrl.level > 0:
                scaling_ctrl.max_feas = max(kkt_result.kkt_resi)
                scaling_ctrl.rel_gap = kkt_result.pd_gap

            if adjust_sigma_yes:
                use_feas_org = use_feas_org or max(kkt_result.kkt_resi) < params.tol_feasorg
                pd_ratio = kkt_result.org_pd_ratio if use_feas_org else kkt_result.pd_ratio
                sigma_new = penalty_manager.get_updated_value(sigma, pd_ratio)
                factor = sigma_new / sigma if sigma != 0 else 1.0
                sigma = sigma_new
                if abs(factor - 1.0) > 1e-12:
                    logging.debug(f"  -> Adjusting Sigma: {sigma:.2e} (Factor: {factor:.2f})")
                    state.alpha /= factor
                    state.beta /= factor
                    state.c_scaled /= factor

        elapsed = time.perf_counter() - total_start_time
        if elapsed >= params.time_limit:
            break
    
    history_manager.end()

    # Final logging
    if is_converged:
        logging.info(f"Converged at iter {it} with KKT error {kkt_result.validation:.2e}\n")
    else:
        if it >= params.maxit:
            logging.info(f"Reached maximum iterations ({params.maxit}) with KKT error {kkt_result.validation:.2e}\n")
        elif elapsed >= params.time_limit:
            logging.info(f"Stopped at iter {it} due to time limit with KKT error {kkt_result.validation:.2e}\n")
        else:
            raise RuntimeError("Solver terminated unexpectedly.")

    # Prepare output
    var_out = copy.deepcopy(var)
    var_out.phi = state.phi
    var_out.q = state.q
    var_out.z = state.z
    var_out.alpha = sigma * state.alpha
    var_out.beta = sigma * state.beta
    var_out.cScale = state.c_scale
    var_out.dScale = state.d_scale
    var_out.D = context.D
    var_out.E = context.E
    sigma_out = sigma / scaling_ctrl.sigma_scale

    return var_out, history_manager, sigma_out

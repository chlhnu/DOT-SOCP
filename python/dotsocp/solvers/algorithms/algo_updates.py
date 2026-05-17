from typing import Optional, Callable
from contextlib import contextmanager
import numexpr as ne

from .utils.admm_utils import RunHistoryManager
from .utils.operators import solve_poisson
from .utils.dataclass import (
    SolverParams,
    ProblemContext,
    IterationState,
)

@contextmanager
def _noop_context():
    yield


def _get_timer(history_manager: Optional[RunHistoryManager]) -> Callable:
    if history_manager is None:
        return lambda **_: _noop_context()
    return history_manager.timer


def inexact_palm_update(
    state: IterationState,
    context: ProblemContext,
    params: SolverParams,
    history_manager: Optional[RunHistoryManager] = None,
):
    # Get timer for subproblem steps
    timer = _get_timer(history_manager)

    # Get context
    scale_d = context.E / state.d_scale
    backend = context.backend
    weight = 1.0 if context.weight is None else context.weight

    # Get parameters
    tau = params.tau

    # Update steps
    with timer(tag="Step phi"):
        # rhs = context.At @ (weight * state.q - state.alpha) + state.c_scaled
        rhs = context.At @ ne.evaluate(
            "weight * q - alpha",
            local_dict={"weight": weight, "q": state.q, "alpha": state.alpha}
        ) + state.c_scaled
        state.phi = solve_poisson(context.kernel, rhs)
        state.tmp_q = context.grad @ state.phi

    with timer(tag="Step z"):
        # backend.proj_soc(state.z, state.z2 - state.beta)
        backend.proj_soc(state.z, ne.evaluate("z2 - beta", local_dict={"z2": state.z2, "beta": state.beta}))

    with timer(tag="Step q"):
        # backend.oper_bfd_conj(state.q2, state.z + state.beta, context.scale_bf)
        backend.oper_bfd_conj(state.q2, ne.evaluate("z + beta", local_dict={"z": state.z, "beta": state.beta}), context.scale_bf)
        # state.q = (weight * (state.tmp_q + state.alpha) + state.q2) * context.diag_q_inv
        ne.evaluate(
            "(weight * (tmp_q + alpha) + q2) * diag_q_inv",
            local_dict={
                "weight": weight,
                "tmp_q": state.tmp_q,
                "alpha": state.alpha,
                "q2": state.q2,
                "diag_q_inv": context.diag_q_inv,
            },
            out=state.q,
        )
        backend.oper_bfd(state.z2, state.q, context.scale_bf, scale_d)

    with timer(tag="Step multipliers"):
        # state.resi_alpha = state.tmp_q - weight * state.q
        # state.alpha += tau * state.resi_alpha
        ne.evaluate("tmp_q - weight * q", local_dict={"tmp_q": state.tmp_q, "weight": weight, "q": state.q}, out=state.resi_alpha)
        ne.evaluate("alpha + tau * resi_alpha", local_dict={"alpha": state.alpha, "tau": tau, "resi_alpha": state.resi_alpha}, out=state.alpha)

        # state.resi_beta = state.z - state.z2
        # state.beta += tau * state.resi_beta
        ne.evaluate("z - z2", local_dict={"z": state.z, "z2": state.z2}, out=state.resi_beta)
        ne.evaluate("beta + tau * resi_beta", local_dict={"beta": state.beta, "tau": tau, "resi_beta": state.resi_beta}, out=state.beta)


def palm_update(
    state: IterationState,
    context: ProblemContext,
    params: SolverParams,
    history_manager: Optional[RunHistoryManager] = None,
):
    # Get timer for subproblem steps
    timer = _get_timer(history_manager)

    # Get context
    scale_d = context.E / state.d_scale
    backend = context.backend
    weight = 1.0 if context.weight is None else context.weight

    # Get parameters
    tau = params.tau

    # Update steps
    with timer(tag="Step q - 1"):
        # backend.oper_bfd_conj(state.q2, state.z + state.beta, context.scale_bf)
        backend.oper_bfd_conj(state.q2, ne.evaluate("z + beta", local_dict={"z": state.z, "beta": state.beta}), context.scale_bf)
        # state.q = (weight * (state.tmp_q + state.alpha) + state.q2) * context.diag_q_inv
        ne.evaluate(
            "(weight * (tmp_q + alpha) + q2) * diag_q_inv",
            local_dict={
                "weight": weight,
                "tmp_q": state.tmp_q,
                "alpha": state.alpha,
                "q2": state.q2,
                "diag_q_inv": context.diag_q_inv,
            },
            out=state.q,
        )
        backend.oper_bfd(state.z2, state.q, context.scale_bf, scale_d)

    with timer(tag="Step phi"):
        # rhs = context.At @ (weight * state.q - state.alpha) + state.c_scaled
        rhs = context.At @ ne.evaluate(
            "weight * q - alpha",
            local_dict={"weight": weight, "q": state.q, "alpha": state.alpha}
        ) + state.c_scaled
        state.phi = solve_poisson(context.kernel, rhs)
        state.tmp_q = context.grad @ state.phi

    with timer(tag="Step z"):
        # backend.proj_soc(state.z, state.z2 - state.beta)
        backend.proj_soc(state.z, ne.evaluate("z2 - beta", local_dict={"z2": state.z2, "beta": state.beta}))

    with timer(tag="Step q - 2"):
        # backend.oper_bfd_conj(state.q2, state.z + state.beta, context.scale_bf)
        backend.oper_bfd_conj(state.q2, ne.evaluate("z + beta", local_dict={"z": state.z, "beta": state.beta}), context.scale_bf)
        # state.q = (weight * (state.tmp_q + state.alpha) + state.q2) * context.diag_q_inv
        ne.evaluate(
            "(weight * (tmp_q + alpha) + q2) * diag_q_inv",
            local_dict={
                "weight": weight,
                "tmp_q": state.tmp_q,
                "alpha": state.alpha,
                "q2": state.q2,
                "diag_q_inv": context.diag_q_inv,
            },
            out=state.q,
        )
        backend.oper_bfd(state.z2, state.q, context.scale_bf, scale_d)

    with timer(tag="Step multipliers"):
        # state.resi_alpha = state.tmp_q - weight * state.q
        # state.alpha += tau * state.resi_alpha
        ne.evaluate("tmp_q - weight * q", local_dict={"tmp_q": state.tmp_q, "weight": weight, "q": state.q}, out=state.resi_alpha)
        ne.evaluate("alpha + tau * resi_alpha", local_dict={"alpha": state.alpha, "tau": tau, "resi_alpha": state.resi_alpha}, out=state.alpha)

        # state.resi_beta = state.z - state.z2
        # state.beta += tau * state.resi_beta
        ne.evaluate("z - z2", local_dict={"z": state.z, "z2": state.z2}, out=state.resi_beta)
        ne.evaluate("beta + tau * resi_beta", local_dict={"beta": state.beta, "tau": tau, "resi_beta": state.resi_beta}, out=state.beta)

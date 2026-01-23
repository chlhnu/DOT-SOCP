from __future__ import annotations

import os
import time
import logging
from dataclasses import dataclass
from typing import List, Optional

import numpy as np

from .algorithms.socp_solver import run_socp_solver
from .algorithms.algo_updates import (
    inexact_palm_update,
    palm_update,
)

from .utils.initialize import initialize_1d, initialize_2d, VariableState, ModelState
from .utils.multilevel import downsample_phi_1d, downsample_phi_2d, downsample_q_2d, jump_next_level
from .utils.scaling import initial_scaling, recover_org_var
from .utils.dot_barrier import ensure_barrier_validity, get_weight_from_barrier

from dotsocp.configs.icons import LOG_ICON


@dataclass
class SolverOutput:
    """Structured output of solver_dotsocp/solve for IDE hints."""

    phi: np.ndarray
    var: VariableState
    model: ModelState
    history: List
    total_time: float

    def last_history(self):
        return self.history[-1] if self.history else None


def solver_dotsocp(rho0_fine, rho1_fine, nt_fine, level_n, opts, method="inpALM", weight=None, barrier_fn=None) -> SolverOutput:
    assert rho0_fine.shape == rho1_fine.shape, "Initial and target densities must have the same shape."
    if (weight is not None) and (barrier_fn is not None):
        raise AssertionError("Only one of weight or barrier_fn can be provided.")

    if barrier_fn is not None:
        rho0_fine, rho1_fine = ensure_barrier_validity(rho0_fine, rho1_fine, barrier_fn)
        weight = get_weight_from_barrier(barrier_fn, nt_fine, rho0_fine.shape[1], rho0_fine.shape[0])

    _resolve_env(weight=weight)

    dim = rho0_fine.squeeze().ndim
    uneven_ot = True if weight is not None else False

    if dim == 1:
        downsample_fn = downsample_phi_1d
        # downsample_q_fn = downsample_q_1d
        initialize_fn = initialize_1d
        assert uneven_ot is False, "Weighted OT not supported for 1D problems."
    elif dim == 2:
        downsample_fn = downsample_phi_2d
        downsample_q_fn = downsample_q_2d
        initialize_fn = initialize_2d
    else:
        raise NotImplementedError(f"Dimension {dim} not supported for multilevel solver.")

    rho0s = [None] * level_n
    rho1s = [None] * level_n
    weights = [None] * level_n
    nts = [None] * level_n
    tols = [None] * level_n

    rho0s[-1] = rho0_fine
    rho1s[-1] = rho1_fine
    weights[-1] = weight.flatten(order='F').reshape(-1, 1) if uneven_ot else None
    nts[-1] = nt_fine

    base_tol = opts.get('tol', 1e-4)
    tols[-1] = base_tol

    tol_factor = -1.0 if base_tol > 0.99e-3 else -0.5
    tol_lower_bound = opts.get('multilevel_tol_lb', None)
    if tol_lower_bound is None:
        tol_lower_bound = 1e-5 if dim == 1 else 1e-4

    for level in range(level_n - 2, -1, -1):
        nts[level] = (nts[level + 1] - 1) // 2 + 1
        tols[level] = max(tols[level + 1] * (2 ** tol_factor), tol_lower_bound)
        rho0s[level] = downsample_fn(rho0s[level + 1])
        rho1s[level] = downsample_fn(rho1s[level + 1])
        rho0s[level] /= np.mean(rho0s[level])
        rho1s[level] /= np.mean(rho1s[level])

        if uneven_ot:
            if barrier_fn is not None:
                # For barrier-type weights, exp(downsampled_log_weight)
                rho0s[level], rho1s[level] = ensure_barrier_validity(rho0s[level], rho1s[level], barrier_fn)
                weights[level] = np.exp(downsample_q_fn(np.log(weights[level + 1]), nts[level + 1], rho0s[level + 1].shape[1], rho0s[level + 1].shape[0]))
            else:
                weights[level] = downsample_q_fn(weights[level + 1], nts[level + 1], rho0s[level + 1].shape[1], rho0s[level + 1].shape[0])

    var: VariableState
    model: ModelState
    time_ml = []
    histories = []

    if method.lower() == "inpalm":
        update_fn = inexact_palm_update
    elif method.lower() == "palm":
        update_fn = palm_update
    else:
        raise NotImplementedError(f"Method {method} not supported.")

    total_start = time.perf_counter()
    elapsed_time = 0.0
    time_limit = opts.get('time_limit', 3600.0)
    time_limit = time_limit if time_limit > 0 else float('inf') # time_limit <= 0 means no time limit
    last_level_kkt = None

    for level in range(level_n):
        current_nt = nts[level]

        if dim == 1:
            current_ny, current_nx = None, rho0s[level].shape[0]
            grid_msg = f"{current_nx}, {current_nt}"
        elif dim == 2:
            current_ny, current_nx = rho0s[level].shape[0], rho0s[level].shape[1]
            grid_msg = f"{current_ny}, {current_nx}, {current_nt}"

        logging.info(
            "{icon}  Entering level {lvl}/{total} | Grid ({grid_msg})\n{sep}".format(
                icon=LOG_ICON['stage'],
                sep='-' * 60,
                lvl=level + 1,
                total=level_n,
                grid_msg=grid_msg,
            )
        )

        if level == 0:
            var, model = initialize_fn(rho0s[level], rho1s[level], current_nt, weight=weights[level])
        else:
            var, model = jump_next_level(var, model, rho0s[level], rho1s[level], current_nt, weight=weights[level])

        scaling_yes = opts.get('scaling', True)
        var, model = initial_scaling(var, model, scaling_yes, last_level_kkt)
        opts_level = opts.copy()
        opts_level['tol'] = tols[level]
        
        remaining_time = time_limit - elapsed_time
        # warning: time_limit <= 0 means no time limit and remaining time <= 0 means time exceeded
        opts_level['time_limit'] = remaining_time if remaining_time > 0 else 1e-6
        
        var, run_history, sigma_final = run_socp_solver(var, model, opts_level, update_fn=update_fn)
        histories.append(run_history)

        if 'sigma' in opts:
            opts['sigma'] = 10 ** (np.log10(opts['sigma'] * sigma_final) / 2)

        if run_history.kkt_errors.size > 0:
            last_level_kkt = run_history.kkt_errors[-1, :]


        var = recover_org_var(var)

        time_ml.append(time.perf_counter() - total_start)
        elapsed_time = time_ml[-1]

    total_time = time.perf_counter() - total_start
    logging.info(
        "{icon}  Solver completed in {time:.2f} sec\n{sep}\n".format(
            icon=LOG_ICON['done'],
            sep='=' * 60,
            time=total_time,
        )
    )

    return SolverOutput(
        phi=var.phi,
        var=var,
        model=model,
        history=histories,
        total_time=total_time,
    )


_AUTO_DOTSOCP_CHECK_RHO_B: Optional[str] = None


def _resolve_env(weight: Optional[np.ndarray] = None) -> None:
    # -- Set default value for env `DOTSOCP_CHECK_RHO_B`
    # A pragmatic heuristic.
    # In the runs with highly concentrated densities and barrier-like weights,
    # this KKT error may look stubborn even when other residuals improve, possibly due to discretization effects on the surface of the barrier.
    # If this is undesired, set the env var explicitly.
    global _AUTO_DOTSOCP_CHECK_RHO_B

    env_key = "DOTSOCP_CHECK_RHO_B"
    current = os.getenv(env_key)

    if current is not None and _AUTO_DOTSOCP_CHECK_RHO_B is None:
        return
    if current is not None and _AUTO_DOTSOCP_CHECK_RHO_B is not None and current != _AUTO_DOTSOCP_CHECK_RHO_B:
        return

    is_barrier_like = weight is not None and np.max(weight) / (1.0 + np.min(weight)) > 1e4
    new_value = "0" if is_barrier_like else "1"
    os.environ[env_key] = new_value
    _AUTO_DOTSOCP_CHECK_RHO_B = new_value

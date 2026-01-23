from __future__ import annotations

import logging
from dataclasses import dataclass
from typing import Callable, Literal

import numpy as np

from dotsocp.solvers import SolverOutput, solver_dotsocp
from dotsocp.visualizers.runhistory import MultilevelRunHistoryManager
from dotsocp.presets.density1d.generator import (
    get_example as get_example_1d,
    get_predefined_examples as list_density_1d,
)
from dotsocp.presets.density2d.generator import (
    get_example as get_example_2d,
    get_predefined_examples as list_density_2d,
)
from dotsocp.presets.weight2d.weight_generator import (
    get_weight as get_weight_2d,
    get_predefined_weights as list_weight_2d,
)
from dotsocp.presets.weight2d.barrier_generator import (
    get_barrier_fn as get_barrier_2d,
    get_predefined_barriers as list_barrier_2d,
)
from dotsocp.visualizers.recovery import (
    recover_density, recover_movement, recover_q,
)
from dotsocp.configs.icons import LOG_ICON


__all__ = [
    # Print presets
    "list_presets",
    "load_presets",
    # Configs
    "GridConfig",
    "SolverConfig",
    "PresetConfig",
    # Solution container
    "DotSolution",
    # Solver
    "solve",
    # Helpers
    # "recover_density",
    # "recover_movement",
    # "recover_q",
]


@dataclass
class GridConfig:
    nt: int
    nx: int
    ny: int | None = None

    def resolved(self) -> tuple[int, int, int]:
        """Return (nt, nx, ny_resolved)."""
        ny_resolved = self.ny if self.ny is not None else self.nx
        return self.nt, self.nx, ny_resolved


@dataclass
class SolverConfig:
    levels: int = 3
    tol: float = 1e-4
    maxit: int = 3000
    time_limit: float = 3600.0
    method: Literal["inpALM", "pALM"] = "inpALM"
    scaling: bool = True
    sigma: float = 1.0
    tau: float = 1.9


@dataclass
class PresetConfig:
    dim: Literal[1, 2]
    density_name: str | None = None
    weight_name: str | None = None
    barrier_name: str | None = None
    lower_bound: float = 0.0

    def __post_init__(self) -> None:
        if self.dim == 1 and (self.weight_name or self.barrier_name):
            raise ValueError("Weight and barrier are only supported for 2D problems.")
        if self.weight_name and self.barrier_name:
            raise ValueError("Specify only one of weight_name or barrier_name.")


@dataclass
class DotSolution:
    """Unified DOT solution container returned by solve()."""
    rho: np.ndarray | None
    Ex: np.ndarray  | None = None
    Ey: np.ndarray  | None = None # 2D only
    q0: np.ndarray  | None = None
    bx: np.ndarray  | None = None
    by: np.ndarray  | None = None # 2D only


def list_presets(kind: Literal["density", "weight", "barrier"], dim: Literal[1, 2]) -> list[str]:
    """Return available preset names."""
    if dim not in (1, 2):
        raise ValueError("Dimension must be 1 or 2.")
    
    if kind == "density":
        return list_density_1d() if dim == 1 else list_density_2d()

    elif kind == "weight":
        if dim == 2:
            return list_weight_2d()
        elif dim == 1:
            raise NotImplementedError("Weight presets for 1D problems are not supported yet.")
    
    elif kind == "barrier":
        if dim == 2:
            return list_barrier_2d()
        elif dim == 1:
            raise NotImplementedError("Barrier presets for 1D problems are not supported yet.")
    
    raise ValueError(f"Unknown preset kind '{kind}'.")


def load_presets(grid: GridConfig, preset: PresetConfig) -> tuple[np.ndarray, np.ndarray, np.ndarray | None, Callable | None]:
    """Load densities and optional weight/barrier according to the preset config."""
    nt, nx, ny = grid.resolved()
    if preset.dim == 1:
        rho0, rho1 = get_example_1d(preset.density_name or "gaussian", nx, lower_bound=preset.lower_bound)
        weight = None
        barrier_fn = None
    else:
        rho0, rho1 = get_example_2d(preset.density_name or "example1", nx, ny, lower_bound=preset.lower_bound)
        weight = get_weight_2d(preset.weight_name, nt, nx, ny) if preset.weight_name else None
        barrier_fn = get_barrier_2d(preset.barrier_name) if preset.barrier_name else None
    return rho0, rho1, weight, barrier_fn


def solve(
    grid: GridConfig,
    solver: SolverConfig,
    preset: PresetConfig,
    logger: logging.Logger | None = None,
) -> tuple[DotSolution, MultilevelRunHistoryManager]:
    """High-level entry for running DOT-SOCP with presets."""
    _ensure_logger(logger)

    rho0, rho1, weight, barrier_fn = load_presets(grid, preset)

    opts = {
        "tol": solver.tol,
        "maxit": solver.maxit,
        "sigma": solver.sigma,
        "tau": solver.tau,
        "scaling": solver.scaling,
        "time_limit": solver.time_limit,
    }

    nt, nx, ny = grid.resolved()
    dim = 1 if preset.dim == 1 else 2
    logging.info(
        "%s  Running DOT-SOCP | dim=%sd | grid=(%s) | levels=%d | density=%s | delta=%s | weight=%s | barrier=%s\n",
        LOG_ICON.get("start", "[run]"),
        dim,
        f"{nx}x{nt}" if dim == 1 else f"{nx}x{ny}x{nt}",
        solver.levels,
        preset.density_name or ("gaussian" if dim == 1 else "example1"),
        preset.lower_bound,
        preset.weight_name or "-",
        preset.barrier_name or "-",
    )

    output: SolverOutput = solver_dotsocp(
        rho0,
        rho1,
        nt,
        solver.levels,
        opts,
        method=solver.method,
        weight=weight,
        barrier_fn=barrier_fn,
    )

    # Recover DOT solution
    rho = recover_density(output)
    movement = recover_movement(output)
    q_vec = recover_q(output)

    if output.model.dim == 1:
        dot_solution = DotSolution(rho=rho, Ex=movement[0], q0=q_vec[0], bx=q_vec[1])
    if output.model.dim == 2:
        dot_solution = DotSolution(rho=rho, Ex=movement[0], Ey=movement[1], q0=q_vec[0], bx=q_vec[1], by=q_vec[2])

    # Multilevel history
    combined_history = MultilevelRunHistoryManager(all_time=output.total_time, histories=output.history)
    # combined_history.print()

    return dot_solution, combined_history


def _ensure_logger(logger: logging.Logger | None) -> None:
    if logger is not None:
        return
    logging.basicConfig(level=logging.INFO, format="%(levelname)s - %(message)s")

import os
import matplotlib.pyplot as plt
from dotsocp import GridConfig, PresetConfig, SolverConfig, load_presets, solve
from dotsocp.visualizers.evolution import show_evolution_1d, show_evolution_2d


def run_case(*,
    dim: int, density: str, weight: str | None = None, barrier: str | None = None, # Problem Setup
    nt: int = 33, nx: int = 129, ny: int | None = None,                            # Grid
    plot_func: str | None = None, plot_cmap: str | None = None,                    # Plot
) -> tuple:
    grid = GridConfig(nt=nt, nx=nx, ny=ny)
    solver = SolverConfig(levels=3, tol=1e-4, maxit=5000, time_limit=1000.0)
    preset = PresetConfig(dim=dim, density_name=density, weight_name=weight, barrier_name=barrier)
    result, history = solve(grid, solver, preset)
    # history.print()

    if dim == 1:
        _ = show_evolution_1d(result.rho, show_func=plot_func or "join", fig_name="Density Evolution")
        plt.show()
    else:
        _, _, _, barrier_fn = load_presets(grid, preset)

        _ = show_evolution_2d(
            result.rho,
            show_func=plot_func or "imshow",
            barrier_fn=barrier_fn,
            fig_name="Density Evolution",
            cmap=plot_cmap,
        )
        plt.show()

    return result, history


def run_dot_1d():
    result, history = run_case(dim=1, density="gaussian", nt=33, nx=129, plot_func="join")
    return result, history

def run_dot_2d(case: str = "ex6"):
    if case == "ex4":
        result, history = run_case(dim=2, density="example4", plot_func="mesh")
    elif case == "ex6":
        os.environ["DOTSOCP_EXAMPLE6_IMG_TYPE"] = "ClassicImages"
        os.environ["DOTSOCP_EXAMPLE6_STITCH1"] = "1,2,3,4"
        os.environ["DOTSOCP_EXAMPLE6_STITCH2"] = "5,6,7,8"
        result, history = run_case(dim=2, density="example6", nt=33, nx=257, plot_func="imshow", plot_cmap="Greys_r")
    else:
        result, history = run_case(dim=2, density="example1", plot_func="imshow", plot_cmap="bone_r")
    return result, history

def run_wdot_weight():
    result, history = run_case(dim=2, density="example1", weight="circular", plot_func="imshow", plot_cmap="bone_r")
    return result, history

def run_wdot_barrier():
    result, history = run_case(dim=2, density="example9", barrier="circle-pillar", plot_func="contourf", plot_cmap="turbo")
    return result, history


if __name__ == "__main__":
    # run_dot_1d()
    run_dot_2d()
    # run_wdot_weight()
    # run_wdot_barrier()

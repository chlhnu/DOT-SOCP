import os
from pathlib import Path
import matplotlib.pyplot as plt
from dotsocp import GridConfig, PresetConfig, SolverConfig, load_presets, solve
from dotsocp.visualizers.evolution import show_evolution_1d, save_evolution_2d

OUTPUT_DIR = Path(__file__).parent / "outputs"


def run_case(*,
    dim: int,density: str, weight: str | None = None, barrier: str | None = None,
    nt: int = 33, nx: int = 129, ny: int | None = None, tol: float = 1e-4,
    plot_func: str | None = None, plot_cmap: str | None = None,
    fig_title: str | None = None, file_name: str | None = None,
) -> tuple:
    grid = GridConfig(nt=nt, nx=nx, ny=ny)
    solver = SolverConfig(levels=3, tol=tol, maxit=5000, time_limit=5000.0)
    preset = PresetConfig(dim=dim, density_name=density, weight_name=weight, barrier_name=barrier)
    result, history = solve(grid, solver, preset)
    # history.print()

    if dim == 1:
        _ = show_evolution_1d(result.rho, show_func=plot_func or "join", fig_name=fig_title or "Density Evolution")
        plt.show()
    else:
        _, _, _, barrier_fn = load_presets(grid, preset)

        _ = save_evolution_2d(
            result.rho, 
            file_name or "evolution-" + "-".join([density, weight or "noweight", barrier or "nobarrier"]) + ".gif",
            show_func=plot_func or "imshow", 
            barrier_fn=barrier_fn,
            fig_name=fig_title or "Density Evolution",
            cmap=plot_cmap
        )

    return result, history


def run_dot_2d_ex4():
    result, history = run_case(
        dim=2, density="example4", plot_func="mesh",
        fig_title="Quartic to Gaussian Mixture", file_name=OUTPUT_DIR / "dot_quartic_to_gaussian.gif",
    )
    return result, history

def run_dot_2d_ex6():
    os.environ["DOTSOCP_EXAMPLE6_IMG_TYPE"] = "ClassicImages"
    os.environ["DOTSOCP_EXAMPLE6_STITCH1"] = "1,2,3,4"
    os.environ["DOTSOCP_EXAMPLE6_STITCH2"] = "5,6,7,8"
    result, history = run_case(
        dim=2, density="example6", plot_func="imshow", plot_cmap="Greys_r",
        nt=33, nx=257,
        fig_title="DOTmark Image Morphing", file_name=OUTPUT_DIR / "dot_dotmark.gif",
    )
    return result, history

def run_dot_2d_ex7():
    result, history = run_case(
        dim=2, density="example7", plot_func="contourf", plot_cmap="turbo",
        nt=129, nx=129,
        fig_title="Gaussian to Dirac", file_name=OUTPUT_DIR / "dot_gaussian_to_dirac.gif",
    )
    return result, history

def run_wdot_weight_circular():
    result, history = run_case(
        dim=2, density="example1", weight="circular", plot_func="imshow", plot_cmap="bone_r",
        nt=129, nx=129,
        fig_title="Weight: Centered Circular", file_name=OUTPUT_DIR / "wdot-circular.gif",
    )
    return result, history

def run_wdot_barrier_circle_pillar():
    result, history = run_case(
        dim=2, density="example9", barrier="circle-pillar", plot_func="contourf", plot_cmap="turbo",
        nt=129, nx=129, tol=1e-5,
        fig_title="Barrier: Circle and Two Pillars", file_name=OUTPUT_DIR / "wdot-circle-pillar.gif",
    )
    return result, history

def run_wdot_barrier_maze():
    result, history = run_case(
        dim=2, density="example9", barrier="maze", plot_func="contourf", plot_cmap="turbo",
        nt=129, nx=129, tol=1e-5,
        fig_title="Barrier: Maze", file_name=OUTPUT_DIR / "wdot-maze.gif",
    )
    return result, history


if __name__ == "__main__":
    run_dot_2d_ex4()
    run_dot_2d_ex6()
    run_dot_2d_ex7()
    run_wdot_weight_circular()
    run_wdot_barrier_circle_pillar()
    run_wdot_barrier_maze()

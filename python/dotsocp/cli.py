from __future__ import annotations

import argparse
import logging

from .interface import (
    GridConfig,
    PresetConfig,
    SolverConfig,
    load_presets,
    list_presets,
    solve,
    DotSolution,
)
from .visualizers.evolution import get_valid_evolution_types, show_evolution_1d, show_evolution_2d


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="DOT-SOCP CLI")
    subparsers = parser.add_subparsers(dest="command", required=True)

    list_parser = subparsers.add_parser("list", help="List available presets")
    list_parser.add_argument("what", choices=["densities", "weights", "barriers"])
    list_parser.add_argument("--dim", type=int, choices=[1, 2], help="Dimension (required for densities)")

    run_parser = subparsers.add_parser("run", help="Run a DOT-SOCP problem from presets")
    run_parser.add_argument("--dim", type=int, choices=[1, 2], default=2, help="Problem dimension")
    run_parser.add_argument("--density", type=str, help="Density preset name")
    run_parser.add_argument("--delta", type=float, default=0.0, help="Add a lower bound of `delta` to the density")
    run_parser.add_argument("--weight", type=str, help="Weight preset name (2D only)")
    run_parser.add_argument("--barrier", type=str, help="Barrier preset name (2D only)")
    run_parser.add_argument("--nt", type=int, default=33, help="Number of time steps")
    run_parser.add_argument("--nx", type=int, default=129, help="Grid size along x")
    run_parser.add_argument("--ny", type=int, default=None, help="Grid size along y (defaults to nx)")
    run_parser.add_argument("--levels", type=int, default=3, help="Number of multilevel refinements")
    run_parser.add_argument("--tol", type=float, default=1e-4, help="Solver tolerance")
    run_parser.add_argument("--maxit", type=int, default=3000, help="Maximum iterations")
    run_parser.add_argument("--time-limit", type=float, default=3600, help="Solver time limit (seconds). None or 0 for no time limit.")
    run_parser.add_argument("--method", choices=["inpALM", "pALM"], default="inpALM", help="Update method")
    run_parser.add_argument("--sigma", type=float, default=1.0, help="Initial sigma")
    run_parser.add_argument("--tau", type=float, default=1.9, help="Initial tau")
    run_parser.add_argument("--no-scaling", action="store_true", help="Disable variable scaling")
    run_parser.add_argument("--plot", action="store_true", help="Show density evolution plot")

    valid_plot_types = get_valid_evolution_types(1) + get_valid_evolution_types(2)
    run_parser.add_argument("--plot-func", type=str, default=None, choices=valid_plot_types, help="Visualization type for density evolution plot")
    run_parser.add_argument("--plot-cmap", type=str, default=None, help="Colormap for density evolution plot. See matplotlib colormaps for valid names (or use 'parula').")

    run_parser.add_argument("--print-history", action="store_true", help="Print running history information of the solver")
    run_parser.add_argument("--plot-history", action="store_true", help="Plot running history curves of the solver")

    return parser


def main(argv: list[str] | None = None) -> None:
    parser = _build_parser()
    args = parser.parse_args(argv)
    logging.basicConfig(level=logging.INFO, format="%(levelname)s - %(message)s")

    if args.command == "list":
        _handle_list(args, parser)
    elif args.command == "run":
        _resolve(args)
        _handle_run(args, parser)


def _resolve(args) -> None:
    # Validate presets
    if args.dim == 1 and (args.weight or args.barrier):
        raise ValueError("Weight and barrier currently are only supported for 2D problems.")
    
    # Validate grid sizes
    if args.levels > 1:
        downsample_factor = 2 ** (args.levels - 1)
        downsample_check = [(x - 1) % downsample_factor == 0 for x in (args.nt, args.nx, args.ny or args.nx) if x is not None]
        if not all(downsample_check):
            raise ValueError(f"For the convenience of {args.levels}-levels multilevel coarsening, nt, nx, ny must satisfy (size - 1) % (2 ** {args.levels - 1}) == 0.")

    # Validate plot-func
    if args.plot and args.plot_func is not None:
        dim = args.dim
        valid_types = get_valid_evolution_types(dim)
        if args.plot_func not in valid_types:
            raise ValueError(f"Invalid plot-func '{args.plot_func}' for dimension {dim}. Must be one of: {', '.join(valid_types)}")
    
    # Validate plot-cmap
    if args.plot and args.plot_cmap is not None:
        if args.plot_cmap.lower() == "parula": # parula colormap defined in MATLAB
            args.plot_cmap = "parula"
            return

        import matplotlib.pyplot as plt
        try:
            plt.get_cmap(args.plot_cmap)
        except ValueError:
            raise ValueError(f"Invalid plot-cmap '{args.plot_cmap}'. See matplotlib colormaps for valid names.")
    
    # Validate tau
    if args.tau < 1e-10 or args.tau >= 2.0:
        raise ValueError("tau must be in the range (0, 2).")


def _handle_list(args, parser) -> None:
    if args.dim not in (1, 2):
        parser.error("--dim must be provided for listing presets (1 or 2)")

    if args.what == "densities":
        names = list_presets("density", dim=args.dim)
    elif args.what == "weights":
        names = list_presets("weight", dim=args.dim)
    else:
        names = list_presets("barrier", dim=args.dim)
    print("\n".join(names))


def _handle_run(args, parser) -> None:
    if args.weight and args.barrier:
        parser.error("Only one of --weight or --barrier can be set.")
    if args.dim == 1 and (args.weight or args.barrier):
        parser.error("Weight/barrier options are only available for 2D problems.")

    grid_cfg = GridConfig(nt=args.nt, nx=args.nx, ny=args.ny)
    solver_cfg = SolverConfig(
        levels=args.levels,
        tol=args.tol,
        maxit=args.maxit,
        time_limit=args.time_limit,
        method=args.method,
        sigma=args.sigma,
        tau=args.tau,
        scaling=not args.no_scaling,
    )
    preset_cfg = PresetConfig(
        dim=args.dim,
        density_name=args.density,
        weight_name=args.weight,
        barrier_name=args.barrier,
        lower_bound=args.delta,
    )

    result: DotSolution
    result, runhistory = solve(grid_cfg, solver_cfg, preset_cfg)

    if args.print_history:
        runhistory.print()
    
    if args.plot_history:
        runhistory.plot()

    if args.plot:
        import matplotlib.pyplot as plt

        if result.rho.ndim == 2:
            _ = show_evolution_1d(result.rho, show_func=args.plot_func, fig_name="Density Evolution")
        else:
            _, _, _, barrier_fn = load_presets(grid_cfg, preset_cfg)
            _ = show_evolution_2d(result.rho, show_func=args.plot_func, barrier_fn=barrier_fn, fig_name="Density Evolution", cmap=args.plot_cmap)
        
        plt.show()


if __name__ == "__main__":
    main()

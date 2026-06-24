DOTSOCP-py
============================================

This folder contains the Python implementation. For the MATLAB version used in the paper's experiments, see ../matlab/README.md.

DOTSOCP-py implements an efficient Second-Order Cone Programming (SOCP) approach for solving Dynamic Optimal Transport (DOT) on staggered grids.

Designed for flexibility and performance, core features include:

- Versatile Problem Solving: Solves DOT (1D/2D) as well as Weighted DOT (2D) problems.
- Fast & Robust Numerics: Uses staggered-grid discretization and multilevel strategies, supporting stable and efficient numerical performance.
- Comprehensive Tooling: Includes CLI, Python API, and visualization modules.

Note: This is research software under active development. Interfaces may change, and feedback is welcome.

Citation
--------
If you use this code in your research, please cite:

- Liang Chen, Youyicun Lin, and Yuxuan Zhou, An efficient second-order cone programming approach for dynamic optimal transport on staggered grid discretization, Mathematical Programming Computation, 2026.

Contact
--------
- E-mail: chl@hnu.edu.cn
- Home page: https://grzy.hnu.edu.cn/site/index/chenliang3

License
--------
DOTSOCP-py is distributed under the GNU Affero General Public License, version 3. See ../LICENSE.

Installation
--------
Ensure you are in the python subdirectory: cd python

Install uv if it is not already installed:

  pip install uv

Create or synchronize the virtual environment and install the CLI tool dotsocp:

  uv sync

Quickstart
--------
Activate the environment:

- Windows (Powershell): .venv/Scripts/activate
- Linux (Bash): source .venv/bin/activate

Solve a transport problem and visualize the result:

- dotsocp run --dim 1 --density gaussian --plot
- dotsocp run --dim 2 --density example1 --plot

List available presets:

- dotsocp list densities --dim 1
- dotsocp list densities --dim 2
- dotsocp list weights --dim 2
- dotsocp list barriers --dim 2

Run weighted DOT or DOT with obstacles:

- Weighted DOT: dotsocp run --dim 2 --density example1 --weight circular --plot
- DOT with Obstacles: dotsocp run --dim 2 --density example9 --barrier maze --plot

Alternatively, invoke the module directly:

  python -m dotsocp <list|run> ...

Examples and Demos
--------
Run lightweight examples with real-time visualization:

  python demos/demos_quick.py

The notebook demos/notebook.ipynb provides a step-by-step guide for 1D/2D DOT and Weighted DOT.
To run it, first install the extra dependencies:

  uv sync --group notebook

Reproduce the GIFs shown above:

  python demos/demos_slow.py

Python API
--------
DOTSOCP-py provides a high-level API to integrate the solver into your own research scripts or notebooks.
The core function solve() returns a DotSolution container and a RunHistoryManager manager.

The basic workflow relies on three configuration objects (GridConfig, SolverConfig, and PresetConfig) passed to the solve() function:

  import matplotlib.pyplot as plt
  import dotsocp
  from dotsocp.visualizers.evolution import show_evolution_2d

  # 1. Define grid size and solver parameters
  grid = dotsocp.GridConfig(nt=33, nx=129, ny=129)
  solver = dotsocp.SolverConfig(levels=3, tol=1e-4, maxit=3000)

  # 2. Select a problem preset
  preset = dotsocp.PresetConfig(dim=2, density_name="example1")

  # 3. Run the solver
  result, history = dotsocp.solve(grid, solver, preset)

  # 4. Analyze and visualize
  history.plot()
  show_evolution_2d(result.rho, show_func="imshow", fig_name="2D Density Evolution")
  plt.show()

---

Thank you for your interest in our work on Dynamic Optimal Transport. We hope this repository serves as a valuable resource for furthering research in this area.

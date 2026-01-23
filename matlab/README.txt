DOTSOCP (MATLAB)
=====================================================

This folder contains the MATLAB implementation. For the Python version, see ../python/README.md.

DOTSOCP implements an efficient Second-Order Cone Programming (SOCP) approach for solving Dynamic Optimal Transport (DOT) on staggered grids.
This folder is used for the numerical experiments in the paper.

Designed for clarity and reproducibility, core features include:

- Self-Contained Solvers: Independent implementations for DOT 1D, DOT 2D, and Weighted DOT 2D, enabling focused usage and easy modification for specific problems.
- Fast, Robust, and Memory-Efficient Numerics: Uses staggered-grid discretization and multilevel strategies for stable and efficient performance, with a low memory footprint.
- Ready-to-Run Experiments: Includes demo scripts, example generators, and lightweight visualization utilities for quick exploration and debugging.

Note
--------------
- This is research software under active development. Interfaces may change, and feedback is welcome

Citation
--------
If you use this code in your research, please cite:

- Liang Chen, Youyicun Lin, and Yuxuan Zhou, An efficient second-order cone programming approach for dynamic optimal transport on staggered grid discretization, arXiv:2505.05424, 2025.

Contact
-------
- E-mail: chl@hnu.edu.cn
- Home page: https://grzy.hnu.edu.cn/site/index/chenliang3

License
-------
DOTSOCP (MATLAB) is distributed under the GNU Affero General Public License, version 3. See ../LICENSE.

Quickstart
----------

DOTSOCP includes three core solvers:

- socp/dot1d: 1-dimensional DOT
- socp/dot2d: 2-dimensional DOT
- socp/wdot2d: 2-dimensional Weighted DOT

These solvers are organized into independent folders to ensure each implementation is self-contained, facilitating focused usage and easy modification for specific problems.

To get started, run the demos:

- demo_dot1d.m: 1-dimensional DOT
- demo_dot2d.m: 2-dimensional DOT
- demo_wdot2d.m: 2-dimensional Weighted DOT

---

Thank you for your interest in our work on Dynamic Optimal Transport. We hope this repository serves as a valuable resource for furthering research in this area.

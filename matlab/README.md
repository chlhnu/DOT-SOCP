# DOTSOCP (MATLAB)

> *Note: This folder contains the **MATLAB implementation**. For the Python version, please see the [../python](../python/) directory.*

[![arXiv](https://img.shields.io/badge/arXiv-2505.05424-b31b1b.svg)](https://arxiv.org/abs/2505.05424)
[![License: AGPL v3](https://img.shields.io/badge/License-AGPL_v3-blue.svg)](../LICENSE)
[![MATLAB 2023b+](https://img.shields.io/badge/MATLAB-2023b+-e16737.svg)](https://www.mathworks.com/products/matlab.html)

This software implements an efficient Second-Order Cone Programming (SOCP) approach for solving Dynamic Optimal Transport (DOT) on staggered grids. 
This folder is used for the numerical experiments in the paper.

Designed for clarity and reproducibility, core features include:

* **Self-Contained Solvers**: Independent implementations for **DOT 1D**, **DOT 2D**, and **Weighted DOT 2D**, enabling focused usage and easy modification for specific problems.
* **Fast, Robust, and Memory-Efficient Numerics**: Uses staggered-grid discretization and multilevel strategies for stable and efficient performance, with a low memory footprint.
* **Ready-to-Run Experiments**: Includes demo scripts, example generators, and lightweight visualization utilities for quick exploration and debugging.

> *Note: This is research software under active development. Interfaces may change, and feedback is welcome.*

## Citation

If you use this code in your research, please cite:

* **Liang Chen, Youyicun Lin, and Yuxuan Zhou**, An efficient second-order cone programming approach for dynamic optimal transport on staggered grid discretization, arXiv:2505.05424, 2025.

## Contact

* E-mail: [chl@hnu.edu.cn](chl@hnu.edu.cn)
* Home page: [https://grzy.hnu.edu.cn/site/index/chenliang3](https://grzy.hnu.edu.cn/site/index/chenliang3)

## Copyright

DOTSOCP (MATLAB) is distributed under the GNU AFFERO GENERAL PUBLIC LICENSE, version 3. A copy can be found in the [LICENSE](../LICENSE).

## Quickstart

**DOTSOCP** includes three core solvers:

* The folder `socp/dot1d` is for 1-dimensional DOT

* The folder `socp/dot2d` is for 2-dimensional DOT

* The folder `socp/wdot2d` is for 2-dimensional Weighted DOT

These solvers are organized into independent folders to ensure each implementation is self-contained, facilitating focused usage and easy modification for specific problems.

To get started, users can run the following simple demos:

* `demo_dot1d.m` is for 1-dimensional DOT

* `demo_dot2d.m` is for 2-dimensional DOT

* `demo_wdot2d.m` is for 2-dimensional Weighted DOT

---

Thank you for your interest in our work on Dynamic Optimal Transport. We hope this repository serves as a valuable resource for furthering research in this area.

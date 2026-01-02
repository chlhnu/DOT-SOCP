%% Demo: solve SOCP for Dynamic Optimal Transport (DOT) 2D

clear;
path(pathdef);
addpath("utils/");
addpath("examples/dot2d/");
addpath("socp/dot2d/");

% KKT error tolerance
tol = 1e-4;

% Discrete grid
nt = 2^5 + 1;
nx = 2^7 + 1;
ny = nx;

% Number of levels
levelN = 3;

%% Set transport problem (Type 1)
% Set <Problem>:
%   "example1"  - Example 5.1
%   "example2"  - Example 5.2
%   "example3"  - Example 5.3
%   "example4"  - Example 5.4
%   "example5"  - Example 5.5
%   "example7"  - Example 5.7
%   "circle"    - Circular densities
%   "DOTmark_4stitch" - 4-Stitch of images from DOTmark
%        Options for "DOTmark_4stitch":
%           "DOTmark_type"       - Type of DOTmark dataset, including "ClassicImages" and "Shapes"
%           "stitch1_indices"    - Indices of images for rho0, list of 4 integers in [1,8], e.g., [1,2,3,4]
%           "stitch2_indices"    - Indices of images for rho1, list of 4 integers in [1,8], e.g., [5,6,7,8]
Problem = "DOTmark_4stitch";

% Set <delta>:
delta = 0;

% Generate densities rho0 and rho1
[rho0, rho1] = get_example(Problem, nx, ny, delta, "DOTmark_type", "ClassicImages", "stitch1_indices", [1,2,3,4], "stitch2_indices", [5,6,7,8]);

%% Set transport problem (Type 2: based on any two images)
% path_to_img1 = "examples/dot2d/resources/DOTmark/ClassicImages/1.png";
% path_to_img2 = "examples/dot2d/resources/DOTmark/ClassicImages/8.png";
% [rho0, rho1] = get_example_from_images(path_to_img1, path_to_img2, nx, ny);

%% Set algorithm
% <algo>:
%   "PALM"          - proximal ALM
%   "inPALM"        - inexact proximal ALM
%   "ALG2"          - ALG2
%   "sGS-inPALM"    - sGS-based inexact proximal ALM
%   "acc-ADMM"      - Accelerated ADMM
%   "acc-sGS-ADMM"  - Accelerated sGS-based ADMM
algo = "inPALM";

%% Solve
opts = {};
opts.tol = tol; % KKT error tolerance
opts.maxit = 3000; % Max iteration

% Solve
[output, timeML, runHistML, runHist] = solver_dotsocp2d(rho0, rho1, nt, levelN, opts, algo);

rho = output.rho;
Ex = output.Ex;
Ey = output.Ey;
q0 = output.q0;
bx = output.bx;
by = output.by;

%% Display
% ---- Evolution ----
show_evolution_2d(rho, "imshow", "Density evolution of " + algo);
% show_evolution_2d(rho, "mesh", "Density evolution of " + algo);
% show_evolution_2d(rho, "contourf", "Density evolution of " + algo);
% show_movement_2d(rho, Ex, Ey, "Density movement of " + algo);

% ---- Violation of rho >= 0 ----
% hist_negative_density(rho, "Histogram: density less than 0 of " + algo);

% ---- Violation of f(q) <= 0 ----
% hist_violation_q_2d(q0, bx, by, "Histogram: f(q) more than 0 of " + algo);

% ---- Mass conservation ----
fprintf(join(repmat("=", 64, 1), "") + "\nCheck mass conservation of %s:\n", algo);
check_massConservation_2d(rho);

% ---- Running history ----
% show_residualCurve(runHistML.kkt, sprintf("KKT errors of %s on all levels", algo), runHistML.kktNames, 'xTime', runHistML.time);
% show_residualCurve(runHist.kkt, sprintf("KKT errors of %s on the last level", algo), runHist.kktNames, 'xIteration', runHist.iter);
% show_residualCurve(runHist.pdGap, sprintf("Primal-dual-gap of %s on the last level", algo), [], 'xIteration', runHist.iter, 'yLabel', 'Relative primal dual gap');

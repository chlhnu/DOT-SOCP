%% Demo: solve SOCP for Weighted Dynamic Optimal Transport (WDOT) 2D

clear;
path(pathdef);
addpath("utils\");
addpath("examples\wdot2d\");
addpath("socp\wdot2d\");

% KKT error tolerance
tol = 1e-5;

% Discrete grid
nt = 2^7 + 1;
nx = 2^7 + 1;
ny = nx;

% Number of levels
levelN = 3;

%% Set Initial and final density
% Set <Problem>:
%   "example1"  - Example 5.1
%   "example2"  - Example 5.2
%   "example3"  - Example 5.3
%   "example4"  - Example 5.4
%   "circle"    - Circular densities
%   "example6"  - Densities of example 5.6
%   "maze14"    - Densities from [Optimal Transport with Proximal Splitting. SIAM Journal on Imaging Sciences, 2014.]

Problem = "example6";
[rho0, rho1] = get_example(Problem, nx, ny);

%% Set weight (Type 1: directly)
% Set <Weight>
%   gene_weight_circle();
%   gene_weight_circleInv();

% weight = gene_weight_circle(nt, nx, ny);

%% Set weight (Type 2: based on barrier)
% Set <barrier>
%   gene_barrier_of_example6(); - Obstacle of Example 5.6
%   gene_barrier_of_maze14();   - Obstacle from [Optimal Transport with Proximal Splitting. SIAM Journal on Imaging Sciences, 2014.]

barrier = gene_barrier_of_example6();
weight = get_weight_by_barrier(nx, ny, nt, barrier);

% Check/Ensure validity of problem
% barrierh = check_barrier_validity(rho0, rho1, barrier);
[rho0, rho1, barrierh] = ensure_barrier_validity(rho0, rho1, barrier);

%% Set algorithm
% Set <algo>:
%   "inPALM"    - Inexact proximal ALM
%   "ALG2"      - ALG2
%   "acc-ADMM"  - Accelerated ADMM
algo = "inPALM";

%% Solver
opts = {};
opts.tol    = tol;
opts.weight = weight;
opts.maxit  = 1e4;

% Solve
if exist("barrier", "var") && ~ isempty(barrier)
    [output, timeML, runHistML, runHist] = solver_wdotsocp2d(rho0, rho1, nt, levelN, opts, algo, barrier);
else
    [output, timeML, runHistML, runHist] = solver_wdotsocp2d(rho0, rho1, nt, levelN, opts, algo);
end

rho = output.rho;
Ex  = output.Ex;
Ey  = output.Ey;
q0  = output.q0;
bx  = output.bx;
by  = output.by;

%% Display
% ---- Evolution ----
if ~ exist("barrierh", "var")
    barrierh = [];
end
show_evolution_2d(rho, "contourf", "Density evolution of " + algo, barrierh);
% show_movement_2d(rho, Ex, Ey, "Density movement of " + algo, barrierh);

% ---- Violation of rho >= 0 ----
% hist_negative_density(rho, "Histogram: density less than 0 of " + algo);

% ---- Violation of f(q) <= 0 ----
% hist_violation_q_2d(q0, bx, by, "Histogram: f(q) more than 0 of " + algo);

% ---- Mass conservation ----
fprintf(join(repmat("=", 64, 1), "") + "\nCheck mass conservation of %s:\n", algo);
check_massConservation_2d(rho);

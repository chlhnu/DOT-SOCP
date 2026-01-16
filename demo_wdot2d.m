%% Demo: Weighted Dynamic Optimal Transport (WDOT) 2D - Quick Start
% This demo solves a weighted 2D optimal transport problem using SOCP.
%
% To see all available presets, run:
%   >> run examples/show_presets_wdot2d
%   >> run examples/show_visual_tools_2d
%
% For detailed examples, see: examples/wdot2d/

clear;
path(pathdef);
addpath("utils/");
addpath("examples/wdot2d/");
addpath("socp/wdot2d/");

% Grid parameters
nt = 33;
nx = 129;
ny = nx;
levelN = 3;

%% Problem Setup
algo = "inpALM";

% -- Type 1 > Generate densities (barrier problems)
problem = "circle-pillar"; % Try: "example8", "example9", "circle-pillar", "maze14"
[rho0, rho1] = get_example(problem, nx, ny);
[weight, barrier] = get_weight_from_barrier(problem, nx, ny, nt);
[rho0, rho1, barrierh] = ensure_barrier_validity(rho0, rho1, barrier);

% -- Type 2 > Generate densities (general weighted problems)
% problem = "example1";   % Try: "example1"-"example4", "circle"
% weight_name = "circle"; % Try: "circle", "circleInv"
% [rho0, rho1] = get_example(problem, nx, ny);
% weight = get_weight(weight_name, nx, ny, nt);
% barrier = []; barrierh = []; % No barrier

%% Solve
opts.tol = 1e-4;        % KKT error tolerance
opts.weight = weight;   % Weight matrix
opts.barrier = barrier; % Barrier function handle
opts.maxit = 1e4;       % Max iteration

[output, timeML, runHistML, runHist] = solver_wdotsocp2d(rho0, rho1, nt, levelN, opts, algo);

%% Visualize Results
show_evolution_2d(output.rho, "contourf", "Density evolution - " + algo, barrierh);
check_massConservation_2d(output.rho);

% For more visualization options, run: 
%   >> run examples/show_visual_tools_2d

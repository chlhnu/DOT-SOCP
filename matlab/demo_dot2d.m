%% Demo: Dynamic Optimal Transport (DOT) 2D - Quick Start
% This demo solves a 2D optimal transport problem using SOCP.
%
% To see all available presets, run:
%   >> run examples/show_presets_dot2d
%   >> run examples/show_visual_tools_2d
%
% For detailed examples, see: examples/dot2d/

clear;
path(pathdef);
addpath("utils/");
addpath("examples/dot2d/");
addpath("socp/dot2d/");

% Grid parameters
nt = 33;
nx = 129;
ny = nx;
levelN = 3;

%% Problem Setup
problem = "example6"; % Try: "example1" - "example7", "circle"
algo = "inpALM";
[rho0, rho1] = get_example(problem, nx, ny);

%% Solve
opts.tol = 1e-4;     % KKT error tolerance
opts.maxit = 3000;   % Max iteration
[output, timeML, runHistML, runHist] = solver_dotsocp2d(rho0, rho1, nt, levelN, opts, algo);

%% Visualize Results
show_evolution_2d(output.rho, "imshow", "Density evolution - " + algo);
check_massConservation_2d(output.rho);

% For more visualization options, run: 
%   >> run examples/show_visual_tools_2d

%% Demo: Dynamic Optimal Transport (DOT) 1D - Quick Start
% This demo solves a 1D optimal transport problem using SOCP.
%
% To see all available presets, run:
%   >> run examples/show_presets_dot1d
%   >> run examples/show_visual_tools_1d
%
% For detailed examples, see: examples/dot1d/

clear;
path(pathdef);
addpath("utils/");
addpath("examples/dot1d/");
addpath("socp/dot1d/");

% Grid parameters
nt = 2^5 + 1;
nx = 2^10 + 1;
levelN = 3;

%% Problem Setup
problem = "gaussian";   % Try: "gaussian", "box"
algo = "inpALM";
[rho0, rho1] = get_example(problem, nx);

%% Solve
opts.tol = 1e-5;
opts.maxit = 3000;
[output, timeML, runHistML, runHist] = solver_dotsocp1d(rho0, rho1, nt, levelN, opts, algo);

%% Display
show_evolution_1d(output.rho, "join");
check_massConservation_1d(output.rho);

% For more visualization options, run: 
%   >> run examples/show_visual_tools_1d

%% Demo: solve SOCP for Dynamic Optimal Transport (DOT) 1D

clear;
path(pathdef);
addpath("utils/");
addpath("examples/dot1d/");
addpath("socp/dot1d/");

% KKT error tolerance
tol = 1e-5;

% Discrete grid
nt = 2^5 + 1;
nx = 2^10 + 1;

% Number of levels
levelN = 3;

% Set transport problem
Problem = "gaussian"; % "gaussian" | "box"

% Set algorithm
algo = "inPALM";

%% Solve
opts = {};
opts.tol = tol;
opts.maxit = 3000; % Maximum number of iterations

% Solve
[rho0, rho1] = get_example(Problem, nx);
[output, timeML, runHistML, runHist] = solver_dotsocp1d(rho0, rho1, nt, levelN, opts, algo);

rho = output.rho;
Ex = output.Ex;
q0 = output.q0;
bx = output.bx;

%% Display
% ---- Evolution ----
show_evolution_1d(rho, "join");
% show_evolution_1d(rho, "tile");

% ---- Violation of rho >= 0 ----
% hist_negative_density(rho, "Histogram: density less than 0 of " + algo);

% ---- Violation of f(q) <= 0 ----
% hist_violation_q_1d(q0, bx, "Histogram: f(q) more than 0 of " + algo);

% ---- Mass conservation ----
fprintf(join(repmat("=", 64, 1), "") + "\nCheck mass conservation of %s:\n", algo);
check_massConservation_1d(rho);

% ---- Running history ----
% show_residualCurve(runHistML.kkt, sprintf("KKT errors of %s on all levels", algo), runHistML.kktNames, 'xTime', runHistML.time);
% show_residualCurve(runHist.kkt, sprintf("KKT errors of %s on the last level", algo), runHist.kktNames, 'xIteration', runHist.iter);
% show_residualCurve(runHist.pdGap, sprintf("Primal-dual-gap of %s on the last level", algo), [], 'xIteration', runHist.iter, 'yLabel', 'Relative primal dual gap');

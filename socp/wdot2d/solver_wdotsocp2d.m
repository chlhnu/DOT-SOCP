function [output, timeML, runHistML, runHist] = solver_wdotsocp2d(rho0, rho1, nt, levelN, opts, method, barrier)
%% Multilevel solver for Weighted DOT-SOCP (2 dimension)
% Input:
%   rho0,   initial density
%   rho1,   terminal density
%   nt,     number of time grid points
%   levelN, number of levels in multilevel scheme
%   opts,   algorithm parameters
%   method, solver type, valid options include:
%       "inPALM" (Inexact Proximal ALM), "ALG2", "acc-ADMM" (Accelerated ADMM)
% 
% Optional input:
%   barrier, function handle of barrier
%
% Output:
%   output,    solution of Weighted DOT, a struct containing fields: rho, Ex, Ey
%   timeML,    timing table for algorithm execution
%   runHistML, algorithm's running history across all levels
%   runHist,   algorithm's running history for the final level
% 
% We will calculate the following KKT errors to evaluate the iterative solution
%   Error of SOCP
%       1 - || A psi - D_w q || / (1 + ||A psi|| + ||q||)
%       2 - || B F q + d - z || / (1 + ||d||)
%       3 - || A^* alpha + c || / (1 + ||c||)
%       4 - || z - Pi_{Q} (z - beta) || / (1 + ||z|| + ||beta||)
%       5 - || F^* B beta + D_w^* alpha || / (1 + ||F^* B beta|| + ||D_w^* alpha||)
%   Error of original DOT
%       1 - || A psi - q || / (1 + ||A psi|| + ||q||)
%       3 - || A^* alpha + c || / (1 + ||c||)
%       6 - || D_{w1}^* alpha_1 - Pi_+ (D_{w1}^* alpha_1 + f(q)) || / (1 + ||D_{w1}^* alpha_1|| + ||f(q)||)
%       7 - || D_{w2}^* alpha_2 - D_{w1}^* g(alpha_1,q) || / (1 + ||D_{w2}^* alpha_2|| + ||D_{w1}^* g(alpha_1,q)||)
% *************************************************************************
% Copyright (c) 2024 by
% Liang Chen, Youyicun Lin, and Yuxuan Zhou
% *************************************************************************

%% Load path
func_full_path = mfilename("fullpath");
[func_dir, ~, ~] = fileparts(func_full_path);

dependent_paths = {'utils', 'algorithms'};

for idx = 1 : length(dependent_paths)
    path = fullfile(func_dir, dependent_paths{idx});
    
    if ~isfolder(path)
        error("path %s not found", path);
    end
    
    addpath(path);
    dependent_paths{idx} = path;
end

cleanup_path = onCleanup(@() rmpath(dependent_paths{:}));

%% Settings
if exist("barrier", "var")
    barrierYes = true;
else
    barrierYes = false;
end

% check levelN
if ~( isnumeric(levelN) && (levelN >= 1) && (levelN == round(levelN)) )
    error("Invalid input at position 4 (Numbers of level in Multilevel strategy)");
end

% Weight
weight = opts.weight;

% check input "method"
if ~exist("method", "var")
    method = "inPALM";
elseif ~ismember(method, ["inPALM", "ALG2", "acc-ADMM"])
    error("Invalid input at position 6 (Solving method)");
end

if (levelN == 1)
    methodName = method + " for Weighted-DOT-SOCP";
else
    methodName = "Multilevel-" + method + " for Weighted-DOT-SOCP";
end

% Whether to check kkt step by step
if ~isfield(opts, "ifCheckStepByStep")
    opts.ifCheckStepByStep = false;
end

% Whether to do sacling
if isfield(opts, "scaling")
    scalingYes = opts.scaling;
else
    scalingYes   = true;
    opts.scaling = true;
end

% Whether to print message
if ~exist("printYes", "var")
    printYes = true;
end

% Max iteration
maxIt = 1e4;

% Stepsize
almStepsize  = 1.9;
alg2Stepsize = 1.0;

% weight
nxs   = cell(levelN, 1);
nys   = cell(levelN, 1);
weights = cell(levelN, 1);
[nys{levelN}, nxs{levelN}] = size(rho0);
weights{levelN} = weight;

% Running history
runHistML = {};
runHistML.kkt   = [];
runHistML.time  = [];
runHistML.iter  = [];
runHistML.pdGap = [];
runHistML.len   = 0;

%% Options
optsML = opts;

% Max iteration
if ~isfield(optsML, "maxit")
    optsML.maxit = maxIt;
end

% Variation of tolerance in Multilevel
if optsML.tol > .99*1e-3
    optsML.tolFactor = -1;
else
    optsML.tolFactor = -.5;
end

tolLowerBound = 1e-4;

% Stepsize
if strcmp(method, "inPALM")
    optsML.tau = almStepsize;
elseif strcmp(method, "ALG2")
    optsML.tau = alg2Stepsize;
end

% Initial sigma
if ~isfield(optsML, 'sigma')
    if isfield(opts, "scaling")
        optsML.sigma = 1;
    else
        optsML.sigma = 0.1;
    end
end

% Time limit
if ~isfield(optsML, 'time_limit')
    optsML.time_limit = 3600;
end

%% Preparation for multilevel
rho0s = cell(levelN, 1);
rho1s = cell(levelN, 1);
nts   = cell(levelN, 1);
tols  = cell(levelN, 1);

rho0s{levelN} = rho0;
rho1s{levelN} = rho1;
nts{levelN}   = nt;

tols{levelN}  = optsML.tol;
tolFactor     = optsML.tolFactor;

for level = levelN-1 : -1 : 1
    nts{level}    = (nts{level+1} + 1) / 2;
    nxs{level}    = (nxs{level+1} + 1) / 2;
    nys{level}    = (nys{level+1} + 1) / 2;
    tols{level}   = max(tols{level+1} * 2^(tolFactor), tolLowerBound);

    % Downsample
    rho0s{level}  = downSample_phi(rho0s{level+1});
    rho1s{level}  = downSample_phi(rho1s{level+1});

    if (barrierYes)
        [rho0s{level}, rho1s{level}] = ensure_barrier_validity(rho0s{level}, rho1s{level}, barrier);
        % weights{level} = get_weight_by_barrier(nxs{level}, nys{level}, nts{level}, barrier);
        weights{level} = downSample_barrier(nts{level+1}, nxs{level+1}, nys{level+1}, weights{level+1});
    else
        weights{level} = downSample_q(nts{level+1}, nxs{level+1}, nys{level+1}, weights{level+1});

        % Normalization
        N = numel(rho0s{level});
        rho0s{level}  = rho0s{level} / (sum(rho0s{level}, 'all') / N);
        rho1s{level}  = rho1s{level} / (sum(rho1s{level}, 'all') / N);
    end
end

%% Multilevel
timeML = cell(levelN + 1, 1);
lastLevelKKT = [];

multilevelClock = tic();

% Initial var and discrete model
[var, model] = initialize(rho0s{1}, rho1s{1}, nts{1});
model.weight = weights{1};

% Multilevel
for level = 1 : levelN
    % Initial scaling
    InitialScaling(var, model, scalingYes, lastLevelKKT);
    
    % Solve
    optsML.tol = tols{level};

    if ismember(method, ["inPALM", "ALG2"])
        [runHist, sigma] = solver_wsocp_inPALM(var, optsML, model);
    elseif strcmp(method, "acc-ADMM")
        [runHist, sigma] = solver_wsocp_accADMM(var, optsML, model);
    end

    % Recover original var
    recoverOrgVar(var);

    % Print message
    if (printYes)
        disp(repmat('=', 1, 64));
        fprintf("Have completed %d-th level (%d, %d, %d)\n", level, nts{level}, size(rho0s{level}));
        disp(var.time);
    end

    % Running history
    timeML{level} = var.time;
    [runHistML, runHist] = catRunHist(runHistML, runHist);

    % Next level
    if level < levelN
        optsML.time_limit = optsML.time_limit - var.time{1,end-1};
        optsML.sigma = 10^(log10(optsML.sigma * sigma) / 2);

        [var, model] = jump_nextLevel(var, model, rho0s{level+1}, rho1s{level+1}, weights{level+1}, nts{level+1});
        lastLevelKKT = runHist.kkt(end, :);
    end
end

multilevelTime = toc(multilevelClock);
disp(repmat('=', 1, 64));
fprintf("Computation time of %s: %.2fs.\n", methodName, multilevelTime);

%% Output
kktLegendNames = {
    '$$\frac{|| \mathcal{A}\psi - D_w q ||}{1 + ||\mathcal{A}\psi|| + ||q||}$$', ...
    '$$\frac{|| \mathcal{B}\mathcal{F} q + d - z ||}{1 + ||d||}$$', ...
    '$$\frac{|| \mathcal{A}^* \alpha + c ||}{1 + ||c||}$$', ...
    '$$\frac{|| z - \Pi_{Q} (z - \beta) ||}{1 + ||z|| + ||\beta||}$$', ...
    '$$\frac{|| \mathcal{F}^* \mathcal{B}^* \beta + D_w^* \alpha ||}{1 + ||\mathcal{F}^* \mathcal{B}^* \beta|| + ||D_w^* \alpha||}$$', ...
    '$$\frac{|| D_{w1}^* \alpha_1 - \Pi_{+} (D_{w1}^* \alpha_1 + f(q)) ||}{1 + ||D_{w1}^* \alpha_1|| + ||f(q)||}$$', ...
    '$$\frac{|| D_{w2}^* \alpha_2 - D_{w1}^* g(\alpha_1,q) ||}{1 + ||D_{w2}^* \alpha_2|| + ||D_{w1}^* g(\alpha_1,q)||}$$'
};

% rho
[rho, Ex, Ey] = recover_RhoE(var, model);

% q = (q0, bx, by)
[q0, bx, by] = recover_q(var, model);

% Solution
output = {};
output.rho = rho;
output.Ex  = Ex;
output.Ey  = Ey;
output.q0  = q0;
output.bx  = bx;
output.by  = by;

% Check mass conservation
tolMass = 1e-2;
conservationYes = check_massConservation(rho, tolMass);
if (~ conservationYes)
    warning("The tolerance of mass conservation constraint is under " + num2str(tolMass));
end

runHistML.method = methodName;
runHistML.kktNames = kktLegendNames;
runHist.method = methodName;
runHist.kktNames = kktLegendNames;

% Running time
timeML{levelN + 1} = record_time(multilevelTime, {'ML_Time'});

end


%% Initial scaling (var/model are handles)
function [] = InitialScaling(var, model, scalingYes, lastLevelKKT)
    h = 1 / numel(var.phi);
    hMean = power(h, 1/3);

    if isempty(lastLevelKKT) || (~isprop(var, "E2"))
        Escale2 = sqrt(2);
    else
        safeguard = 4;
        Escale2 = var.E2 * min(safeguard, max(1 / safeguard, sqrt(lastLevelKKT(1) / lastLevelKKT(2))));
    end

    if scalingYes
        norm_c      = normL2(model.c, h) * sqrt(model.nt);
        norm_d      = sqrt(2);

        adjust      = power(10, mean(log10(model.weight + 1e-10)));
        D           = sqrt(2) * sqrt(hMean) * adjust;
        E           = D / Escale2;
        cScale      = max(1, norm_c * sqrt(hMean) / adjust);
        dScale      = E * norm_d * sqrt(adjust);

        model.c     = model.c    / cScale;
        model.normc = norm_c     / cScale;
        model.normd = norm_d * E / dScale;
        model.grad  = D .* model.grad;
        var.phi     = (1 / dScale)     .* var.phi;
        var.q       = (D / dScale)     .* var.q;
        var.z       = (E / dScale)     .* var.z;
        var.alpha   = (1 / cScale / D) .* var.alpha;
        var.beta    = (1 / cScale / E) .* var.beta;
    else
        cScale      = 1;
        dScale      = 1;
        D           = 1;
        E           = 1;
        model.normc = normL2(model.c, h);
        model.normd = sqrt(2);
    end

    % Scaling factor
    var.cScale = cScale;
    var.dScale = dScale;
    var.D      = D;
    var.E      = E;
    var.E2     = Escale2;
end

%% Recover orignal variables
function recoverOrgVar(var)
    cScale    = var.cScale;
    dScale    = var.dScale;
    D         = var.D;
    E         = var.E;

    var.phi   =    dScale    .* var.phi;
    var.z     = (dScale / E) .* var.z;
    var.q     = (dScale / D) .* var.q;
    var.alpha = (cScale * D) .* var.alpha;
    var.beta  = (cScale * E) .* var.beta;
end

%% Cat running history
function [runHistML, runHist] = catRunHist(runHistML, runHist)
    runHistML.kkt   = cat(1, runHistML.kkt, runHist.kkt);
    runHistML.pdGap = cat(1, runHistML.pdGap, runHist.pdGap);

    if isempty(runHistML.time)
        runHistML.time = runHist.time;
    else
        runHist.time    = runHistML.time(end) + runHist.time;
        runHistML.time  = cat(1, runHistML.time, runHist.time);
    end

    if isempty(runHistML.iter)
        runHistML.iter = runHist.iter;
    else
        runHistML.iter = cat(1, runHistML.iter, runHistML.iter(end) + runHist.iter);
    end

    runHistML.len = runHistML.len + runHist.len;
end

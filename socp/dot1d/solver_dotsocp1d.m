function [output, timeML, runHistML, runHist] = solver_dotsocp1d(rho0, rho1, nt, levelN, opts, method)
%% Multilevel solver for DOT-SOCP (1 dimension)
% Input:
%   rho0,   initial density
%   rho1,   terminal density
%   nt,     number of time grid points
%   levelN, number of levels in multilevel scheme
%   opts,   algorithm parameters
%   method, solver type, valid options include:
%       "inPALM" (Inexact Proximal ALM), "ALG2"
% 
% Output:
%   output,    solution of DOT, a struct containing fields: rho, Ex, q0, bx
%   timeML,    timing table for algorithm execution
%   runHistML, algorithm's running history across all levels
%   runHist,   algorithm's running history for the final level
% 
% We will calculate the following KKT errors to evaluate the iterative solution
%   Error of SOCP
%       1 - || A psi - q || / (1 + ||A psi|| + ||q||)
%       2 - || B F q + d - z || / (1 + ||d||)
%       3 - || A^* alpha + c || / (1 + ||c||)
%       4 - || z - Pi_{Q} (z - beta) || / (1 + ||z|| + ||beta||)
%       5 - || F^* B beta + alpha || / (1 + ||F^* B beta|| + ||alpha||)
%   Error of original DOT
%       1 - || A psi - q || / (1 + ||A psi|| + ||q||)
%       3 - || A^* alpha + c || / (1 + ||c||)
%       6 - || alpha_1 - Pi_+ (f(q) + alpha_1) || / (1 + ||alpha_1|| + ||f(q)||)
%       7 - || alpha_2 - g(alpha_1, q) || / (1 + ||alpha_2|| + ||g(alpha_1, q)||)
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
% check levelN
if ~( isnumeric(levelN) && (levelN >= 1) && (levelN == round(levelN)) )
    error("Invalid input at position 4 (Number of levels in multilevel strategy)");
end

% check input "method"
if ~exist("method", "var")
    method = "inPALM";
elseif ~ismember(method, ["inPALM", "ALG2"])
    error("Invalid input at position 6 (Solving method)");
end

if (levelN == 1)
    methodName = method + " for DOT-SOCP";
else
    methodName = "Multilevel-" + method + " for DOT-SOCP";
end

% Whether to check KKT step by step
if ~isfield(opts, "ifCheckStepByStep")
    opts.ifCheckStepByStep = false;
end

% Whether to do scaling
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
admmMaxIt = 3000;

% Stepsize
almStepsize  = 1.9;
alg2Stepsize = 1.0;

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
if ~isfield(opts, "maxit")
    optsML.maxit = admmMaxIt;
end

% Variation of tolerance in multilevel method
if optsML.tol > .99*1e-3
    optsML.tolFactor = -1;
else
    optsML.tolFactor = -.5;
end

tolLowerBound = 1e-5;

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
    nts{level}    = (nts{level+1} - 1) / 2 + 1;
    tols{level}   = max(tols{level+1} * 2^(tolFactor), tolLowerBound);

    % Downsample
    rho0s{level}  = downSample_phi(rho0s{level+1});
    rho1s{level}  = downSample_phi(rho1s{level+1});

    % Normalization
    N = numel(rho0s{level});
    rho0s{level}  = rho0s{level} / (sum(rho0s{level}, 'all') / N);
    rho1s{level}  = rho1s{level} / (sum(rho1s{level}, 'all') / N);
end

%% Multilevel
timeML = cell(levelN + 1, 1);
lastLevelKKT = [];

multilevelClock = tic();

% Initial var and discrete model
[var, model] = initialize(rho0s{1}, rho1s{1}, nts{1});

% Multilevel
for level = 1 : levelN
    % Initial scaling
    InitialScaling(var, model, scalingYes, lastLevelKKT);
    
    % Solve
    optsML.tol = tols{level};
    [runHist, sigma] = solver_socp_inPALM(var, optsML, model);

    % Recover original var
    recoverOrgVar(var);

    % Print message
    if (printYes)
        disp(repmat('=', 1, 64));
        fprintf("Completed %d-th level (%d, %d)\n", level, nts{level}, size(rho0s{level}, 1));
        disp(var.time);
    end

    % Running history
    timeML{level} = var.time;
    [runHistML, runHist] = catRunHist(runHistML, runHist);

    % Next level
    if level < levelN
        optsML.time_limit = optsML.time_limit - var.time{1,end-1};
        optsML.sigma = 10^(log10(optsML.sigma * sigma) / 2);

        [var, model] = jump_nextLevel(var, model, rho0s{level+1}, rho1s{level+1}, nts{level+1});
        lastLevelKKT = runHist.kkt(end, :);
    end
end

multilevelTime = toc(multilevelClock);
disp(repmat('=', 1, 64));
fprintf("Computation time of %s: %.2fs.\n", methodName, multilevelTime);

%% Output
kktLegendNames = {
    '$$\frac{|| \mathcal{A}\psi - q ||}{1 + ||\mathcal{A}\psi|| + ||q||}$$', ...
    '$$\frac{|| \mathcal{B}\mathcal{F} q + d - z ||}{1 + ||d||}$$', ...
    '$$\frac{|| \mathcal{A}^* \alpha + c ||}{1 + ||c||}$$', ...
    '$$\frac{|| z - \Pi_{Q} (z - \beta) ||}{1 + ||z|| + ||\beta||}$$', ...
    '$$\frac{|| \mathcal{F}^* \mathcal{B}^* \beta + \alpha ||}{1 + ||\mathcal{F}^* \mathcal{B}^* \beta|| + ||\alpha||}$$', ...
    '$$\frac{|| \alpha_1 - \Pi_{+} (\alpha_1 + f(q)) ||}{1 + ||\alpha_1|| + ||f(q)||}$$', ...
    '$$\frac{|| \alpha_2 - g(\alpha_1,q) ||}{1 + ||\alpha_2|| + ||g(\alpha_1,q)||}$$'
};

% rho
[rho, Ex] = recover_RhoE(var, model);

% q = (q0, bx)
[q0, bx] = recover_q(var, model);

% Solution
output = {};
output.rho = rho;
output.Ex  = Ex;
output.q0  = q0;
output.bx  = bx;

% Check mass conservation
tolMass = 1e-2;
conservationYes = check_massConservation(rho, tolMass);
if (~ conservationYes)
    warning("The mass conservation constraint violation exceeds " + num2str(tolMass));
end

% Running history of all levels
runHistML.method = methodName;
runHistML.kktNames = kktLegendNames;

% Running history of the final level
runHist.method = methodName;
runHist.kktNames = kktLegendNames;

% Running time
timeML{levelN + 1} = record_time(multilevelTime, {'ML_Time'});

end


%% Initial scaling
function [] = InitialScaling(var, model, scalingYes, lastLevelKKT)
    h = 1 / numel(var.phi);
    hMean = power(h, 1/2);

    if isempty(lastLevelKKT) || (~isprop(var, "E2"))
        Escale2 = sqrt(2);
    else
        ratio      = sqrt(lastLevelKKT(1) / lastLevelKKT(2));
        lowerRatio = 0.8333;
        if ratio < lowerRatio
            Escale2 = var.E2 * max(1/sqrt(2), ratio / lowerRatio);
        else
            Escale2 = var.E2 * min(sqrt(2), max(1, ratio));
        end
    end

    % Scaling
    if scalingYes
        norm_c      = normL2(model.c, h) * sqrt(model.nt);
        norm_d      = sqrt(2);

        % Scaling factor
        D           = sqrt(2) * sqrt(hMean);
        E           = D / Escale2;
        cScale      = max(1, norm_c * sqrt(hMean) );
        dScale      = E * norm_d;

        % Do scaling
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

%% Recover original variables
function recoverOrgVar(var)
    % Fetch scaling factor
    cScale    = var.cScale;
    dScale    = var.dScale;
    D         = var.D;
    E         = var.E;

    % Recover
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

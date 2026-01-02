function [runHist, sigma] = solver_socp_inPALM(var, opts, model)
%% An inexact proximal ALM for solving the SOCP reformulation of Dynamic Optimal Transport (1 dimension):
%       min <c, \phi> + \delta_{Q}(z)
%       s.t.    A \phi - q = 0,
%               z - B F q  = d.
% *************************************************************************
% Copyright (c) 2024 by
% Liang Chen, Youyicun Lin, and Yuxuan Zhou
% *************************************************************************

%% Params
if ~exist("printYes", "var")
    printYes = false;
end

if ~exist("printDotYes", "var")
    printDotYes = false;
end

if isfield(opts, "checkPrimDualFeas")
    checkPrimDualFeasYes = opts.checkPrimDualFeas;
else
    checkPrimDualFeasYes = true;
end

if isfield(opts, "time_limit")
    time_limit = opts.time_limit;
else
    time_limit = 3600;
end

% Iteration parameters
tau         = opts.tau;
sigma       = opts.sigma;
maxit       = opts.maxit;
tol         = opts.tol;
checkSByS   = opts.ifCheckStepByStep;
lastSigmaIt = -Inf;
updateRule  = [
    1.1, 1.10;
    1.2, 1.15;
    1.5, 1.20;
      2, 1.26;
    2.5, 1.28;
   3.33, 1.32;
      5, 1.35;
     10, 1.40;
     20, 1.60;
     40, 1.80;
     50, 2.00;
];

% Scaling
cScale  = var.cScale;
dScale  = var.dScale;
D       = var.D;
E       = var.E;
scaleBF = E / D;
scaleD  = E / dScale;
use_feasOrg = 0;
tol_feasOrg = 5 * tol;

% Rescaling
if isfield(opts, "scaling") && opts.scaling
    rescale = 1;
else
    rescale = 0;
end

if (rescale == 1)
    firstScaleIter    = 10;
    SecondScaleIter   = 50;
    checkRescaleIters = 100;
    ratioThreshold    = 1.2;
    maxFeas           = Inf;
    relGap            = Inf;
end

%% Initialization
% discrete model
nx      = model.nx;
nt      = model.nt;
h       = 1 / (nx*nt);
A       = model.grad; % Grad
c       = model.c;

% iterative variable
phi     = var.phi;   var.phi   = [];
q       = var.q;     var.q     = [];
z       = var.z;     var.z     = [];
alpha   = var.alpha; var.alpha = [];
beta    = var.beta;  var.beta  = [];

% preprocessing
kernel   = D^2 * initialize_FFTkernel(nt, nx);
diagQInv = 1 ./ oper_q(nx, nt, D, E);

%% Iteration
norm_c = model.normc;
norm_d = model.normd;
alpha  = alpha / sigma;
beta   = beta  / sigma;
c      = c     / sigma;
sigmaScale = 1;

% Running history
kktConst = 1;
runHist  = {};
runHistItems  = 0;
runHist.kkt   = Inf(maxit, 7);
runHist.time  = Inf(maxit, 1);
runHist.iter  = Inf(maxit, 1);
runHist.pdGap = Inf(maxit, 1);

% Stop condition
if (checkPrimDualFeasYes)
    stopCondition = [1,3,6,7];
else
    stopCondition = [1,3,6];
end

% Time
time_lineq      = 0;
time_proj       = 0;
time_q          = 0;
time_multiplier = 0;
time_kkt        = 0;

% Preallocation
z2 = zeros(size(z));
q2 = zeros(size(q));
mexBFd1d(z2, q, nt, nx, scaleBF, scaleD);

clock_total = tic();
for it = 1 : maxit
    % Rescaling
    scaleYes = 0;
    if (rescale >= 3) && (mod(it, checkRescaleIters) == 0)
        [normPhi, normQ, normZ, normAlpha, normBeta] = deal( ...
                normL2(phi, h), normL2(q, h), FnormL2(z, h), sigma * normL2(alpha, h), sigma * FnormL2(beta, h));

        normPhis = max([normPhi, normQ, normZ]);
        normAlps = max([normAlpha, normBeta]);
        ratio = max(normAlps, normPhis) / min(normAlps, normPhis);
        if (ratio > ratioThreshold)
            scaleYes = 1;
        end
    end

    if (rescale == 1) && (maxFeas < 2e-2) && (it >= firstScaleIter) && (relGap < 5e-2) ... % The 1st rescale
            || (rescale == 2) && (maxFeas < 5e-3) && (it >= SecondScaleIter) && (relGap < 1e-2) ... % The 2nd rescale
            || scaleYes % The 3rd, 4-th, ..., n-th rescale

        if (~ scaleYes)
            [normPhi, normQ, normZ, normAlpha, normBeta] = deal( ...
                normL2(phi, h), normL2(q, h), FnormL2(z, h), sigma * normL2(alpha, h), sigma * FnormL2(beta, h));

            normPhis = max([normPhi, normQ, normZ]);
            normAlps = max([normAlpha, normBeta]);
        end

        % Scaling
        %   Rescaling factor
        dScale2 = normPhis;
        cScale2 = normAlps;

        %   Scale parameters
        sigma   = sigma  * (cScale2 / dScale2);
        c       = c      * dScale2 / cScale2^2;
        norm_c  = norm_c / cScale2;
        norm_d  = norm_d / dScale2;

        %   scale var
        alpha   = alpha  * dScale2 / cScale2^2;
        beta    = beta   * dScale2 / cScale2^2;
        q       = q      / dScale2;
        z       = z      / dScale2;

        %   record scaling factor
        dScale  = dScale2 * dScale;
        cScale  = cScale2 * cScale;
        scaleD  = E / dScale;
        sigmaScale = sigmaScale * (cScale2 / dScale2);

        %   update temp var
        mexBFd1d(z2, q, nt, nx, scaleBF, scaleD);

        rescale = rescale + 1;
    end

    % step phi
    clock_lineq = tic();
    phi = oper_poisson(kernel, reshape(A' * (q - alpha) + c, nx, nt));
    time_lineq  = time_lineq + toc(clock_lineq);
    
    % step z
    clock_proj = tic();
    mexProjSoc(z, z2 - beta);
    time_proj = time_proj + toc(clock_proj);
    
    % step q
    clock_q = tic();
    tmp_q = A * phi;
    mexBFdConj1d(q2, z + beta, nt, nx, scaleBF);
    q = (tmp_q + alpha + q2) .* diagQInv;
    time_q = time_q + toc(clock_q);
    
    % step alpha, beta
    clock_multiplier = tic();
    resi_alpha = tmp_q - q;
    mexBFd1d(z2, q, nt, nx, scaleBF, scaleD);
    resi_beta  = z - z2;
    alpha      = alpha + tau * resi_alpha;
    beta       = beta  + tau * resi_beta;
    time_multiplier = time_multiplier + toc(clock_multiplier);
    
    % kkt
    clock_kkt = tic();
    adjustSigmaYes = IfAdjustSigma(it, lastSigmaIt);
    check_kkt_yes = checkSByS || adjustSigmaYes || (it == maxit) || (toc(clock_total) > time_limit);
    if check_kkt_yes
        % Precomputation
        %   temp
        mexBFdConj1d(q2, beta, nt, nx, scaleBF);
        %   norm
        norm_q      = normL2(q, h);
        norm_z      = FnormL2(z, h);
        norm_Aphi   = normL2(tmp_q, h);
        norm_alpha  = sigma * normL2(alpha, h);
        norm_beta   = sigma * FnormL2(beta, h);
        norm_FBbeta = sigma * normL2(q2, h);

        % KKT residuals
        primFea1   = normL2(resi_alpha, h);
        primFea2   = FnormL2(resi_beta, h);
        dualFea1   = sigma * normL2(A'*alpha - c, h);
        dualFea2   = sigma * normL2(q2 + alpha, h);

        mexProjSoc(z2, z - sigma * beta);
        complem    = FnormL2(z - z2, h);
        mexBFd1d(z2, q, nt, nx, scaleBF, scaleD);
        
        [dotcomplem, normRho, norm_rhoFq, mRhoB, normM, normRhoB] = compute_kkt_dot_complement(q, alpha, z2, sigma, h, nt, nx, var.qInd, cScale, dScale, D, E);
        
        % Relative KKT residuals
        KKTResiOrg = [
            primFea1   / (kktConst * D / dScale + norm_Aphi + norm_q), ...
            primFea2   / (kktConst * E / dScale + norm_d), ...
            dualFea1   / (kktConst / cScale + norm_c), ...
            complem    / (kktConst * E / dScale + norm_z + norm_beta), ...
            dualFea2   / (kktConst / cScale / D + norm_FBbeta + norm_alpha), ...
            dotcomplem / (kktConst + normRho + norm_rhoFq), ...
            mRhoB      / (kktConst + normM + normRhoB)
        ];
        KKTResi = [
            primFea1   / (kktConst + norm_Aphi + norm_q), ...
            primFea2   / (kktConst + norm_d), ...
            dualFea1   / (kktConst + norm_c), ...
            complem    / (kktConst + norm_z + norm_beta), ...
            dualFea2   / (kktConst + norm_FBbeta + norm_alpha)
        ];

        % Primal-Dual gap
        priVal   = (sigma * cScale * dScale * h) * dot(q, alpha);
        dualVal  = (sigma * cScale * dScale * h) * dot(c, phi);
        pdGap    = abs(priVal - dualVal) / (1 + abs(priVal) + abs(dualVal));

        % Running history
        runHistItems = runHistItems + 1;
        runHist.kkt(runHistItems, :) = KKTResiOrg;
        runHist.time(runHistItems)   = toc(clock_total);
        runHist.iter(runHistItems)   = it;
        runHist.pdGap(runHistItems)  = pdGap;

        % Print message
        if (printYes)
            fprintf("PrimVal: %.2E. DualVal: %.2E. Prim-Dual gap: %.2E. PrimFeas: %.2E. DualFeas: %.2E. Complement: %.2E.\n", ...
                priVal, dualVal, pdGap, max(KKTResiOrg([1,2])), max(KKTResiOrg([3,5])), KKTResiOrg(4));
        end
        if (printDotYes)
            fprintf("PrimVal: %.2E. DualVal: %.2E. Prim-Dual gap: %.2E. DualFeas: %.2E. PrimFeas: %.2E. Compleme: %.2E. PrimDualFeas: %.2E.\n", ...
                priVal, dualVal, pdGap, KKTResiOrg(1), KKTResiOrg(3), KKTResiOrg(6), KKTResiOrg(7));
        end

        % Stop criterion
        if max(KKTResiOrg(stopCondition)) < tol || ...
                toc(clock_total) > time_limit
            break;
        end

        % Use original KKT
        if max(KKTResi) < tol_feasOrg
            use_feasOrg = 1;
        end
    
        % Update Lagrangian parameter
        if adjustSigmaYes
            lastSigmaIt = it;
            
            if use_feasOrg
                resiPri  = max(KKTResiOrg([1,2]));
                resiDual = max(KKTResiOrg([3,5]));
            else
                resiPri  = max(KKTResi([1,2]));
                resiDual = max(KKTResi([3,5]));
            end

            [sigma, factor] = adjust_lagrangianParam(sigma, resiPri / resiDual, updateRule);

            if factor ~= 1
                alpha  = alpha / factor;
                beta   = beta  / factor;
                c      = c     / factor;
            end
        end

        % Information for rescaling
        if (rescale > 0)
            maxFeas = max(KKTResi);
            relGap  = pdGap;
        end
    end
    time_kkt = time_kkt + toc(clock_kkt);
end
time_total = toc(clock_total);

%% output
var.name = 'Inexact Proximal ALM';

% Iterative var
var.phi = phi;
var.q = q;
var.z = z;
var.alpha = sigma * alpha;
var.beta = sigma * beta;

% Time
times = [time_lineq, time_proj, time_q, time_multiplier, time_kkt, time_total, it];
names = {'Step_1_1_FFT', 'Step_1_2_ProjSOC', 'Step_2_Q_Step', 'Step_3_Multiplier', 'KKT', 'Total_Time', 'Iters'};
var.time = record_time(times, names);

% Running history
runHist.len = runHistItems;
runHist.kkt(runHistItems+1 : end, :)   = [];
runHist.time(runHistItems+1 : end)     = [];
runHist.iter(runHistItems+1 : end)     = [];
runHist.pdGap(runHistItems+1 : end)    = [];

% Scaling factor
var.cScale = cScale;
var.dScale = dScale;
var.D      = D;
var.E      = E;

% Recover sigma
sigma = sigma / sigmaScale;

end

function flag = IfAdjustSigma(iter, last_iter)
    %% Whether to adjust Lagrangian param
    passedIters = iter - last_iter;
    flag = 0;
    
    if iter < 20 && passedIters >= 3
        flag = 1;
    elseif iter < 50 && passedIters >= 6
        flag = 1;
    elseif iter < 100 && passedIters >= 10
        flag = 1;
    elseif iter < 200 && passedIters >= 15
        flag = 1;
    elseif iter < 500 && passedIters >= 25
        flag = 1;
    elseif passedIters >= 40
        flag = 1;
    end
end

function [var_out, runHist, sigma] = solver_socp_inPALM(var, opts, model)
%% An inexact proximal ALM for solving the SOCP reformulation of Dynamic Optimal Transport:
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
ny      = model.ny;
nt      = model.nt;
h       = 1 / (nx*ny*nt);
A       = model.grad; % Grad
At      = transpose(A); % -Div
c       = model.c;

% iterative variable
phi     = var.phi;
q       = var.q;
z       = var.z;
alpha   = var.alpha;
beta    = var.beta;

% preprocessing
kernel   = D^2 * initialize_FFTkernel(nt, nx, ny);
diagQInv = 1 ./ oper_q(ny, nx, nt, D, E);

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
z2 = zeros((nt-1)*nx*ny, 10);
q2 = zeros((nt-1)*nx*ny + nt*(nx-1)*ny + nt*nx*(ny-1), 1);
projZBta = zeros((nt-1)*nx*ny, 10);
mexBFd(z2, q, nt, nx, ny, scaleBF, scaleD);

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
        mexBFd(z2, q, nt, nx, ny, scaleBF, scaleD);

        rescale = rescale + 1;
    end

    % step phi
    clock_lineq = tic();
    phi = oper_poisson3dim(kernel, reshape(At * (q - alpha) + c, ny, nx, nt));
    time_lineq  = time_lineq + toc(clock_lineq);
    
    % step z
    clock_proj = tic();
    mexProjSoc(z, z2 - beta);
    time_proj = time_proj + toc(clock_proj);
    
    % step q
    clock_q = tic();
    tmp_q = A * phi;
    mexBFdConj(q2, z + beta, nt, nx, ny, scaleBF);
    q = (tmp_q + alpha + q2) .* diagQInv;
    time_q = time_q + toc(clock_q);
    
    % step alpha, beta
    clock_multiplier = tic();
    resi_alpha = tmp_q - q;
    mexBFd(z2, q, nt, nx, ny, scaleBF, scaleD);
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
        mexBFdConj(q2, beta, nt, nx, ny, scaleBF);
        mexProjSoc(projZBta, z - sigma * beta);
        %   norm
        norm_q      = normL2(q, h);
        norm_z      = FnormL2(z, h);
        norm_Aphi   = normL2(tmp_q, h);
        norm_alpha  = sigma * normL2(alpha, h);
        norm_beta   = sigma * FnormL2(beta, h);
        norm_FBbeta = sigma * normL2(q2, h);
        %   alpha1
        qInd  = var.qInd;
        rhoT = (sigma * cScale * D) * alpha(1 : qInd.bx-1);
        rhoFq = rhoT + (dScale / D) * q(1 : (nt-1)*ny*nx) + sum(((dScale / E) * z2(:, 2:9)).^2, 2) / 4;
        rhoFq(rhoFq < 0) = 0;
        normRho = normL2(rhoT, h);
        norm_rhoFq  = normL2(rhoFq, h);
        %   alpha2
        rho = movmean(cat(3, zeros(ny, nx), reshape(rhoT, ny, nx, nt-1), zeros(ny, nx)), 2, 3, "Endpoints", "discard");
        rhoBx = (dScale / D) * ( reshape(movmean(rho, 2, 2, "Endpoints", "discard"), [], 1) .* q(qInd.bx : qInd.by-1) );
        rhoBy = (dScale / D) * ( reshape(movmean(rho, 2, 1, "Endpoints", "discard"), [], 1) .* q(qInd.by : end) );
        mx = (sigma * cScale * D) * alpha(qInd.bx : qInd.by-1);
        my = (sigma * cScale * D) * alpha(qInd.by : end);
        normM = sqrt(normL2(mx, h)^2 + normL2(my, h)^2);
        normRhoB = sqrt(normL2(rhoBx, h)^2 + normL2(rhoBy, h)^2);

        % KKT residuals
        primFea1   = normL2(resi_alpha, h);
        primFea2   = FnormL2(resi_beta, h);
        dualFea1   = sigma * normL2(At*alpha - c, h);
        complem    = FnormL2(z - projZBta, h);
        dualFea2   = sigma * normL2(q2 + alpha, h);
        dotcomplem = normL2(rhoT - rhoFq, h);
        mRhoB      = sqrt(normL2(mx - rhoBx, h)^2 + normL2(my - rhoBy, h)^2);
        
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
var_out = var;
var_out.name = 'Inexact Proximal ALM';

% Iterative var
var_out.phi = phi;
var_out.q = q;
var_out.z = z;
var_out.alpha = sigma * alpha;
var_out.beta = sigma * beta;

% Time
times = [time_lineq, time_proj, time_q, time_multiplier, time_kkt, time_total, it];
names = {'Step1-1 FFT', 'Step1-2 ProjSOC', 'Step2 q', 'Step3 multiplier', 'KKT', 'Total time', 'Iters'};
var_out.time = record_time(times, names);

% Running history
runHist.len = runHistItems;
runHist.kkt(runHistItems+1 : end, :)   = [];
runHist.time(runHistItems+1 : end)     = [];
runHist.iter(runHistItems+1 : end)     = [];
runHist.pdGap(runHistItems+1 : end)    = [];

% Scaling factor
var_out.cScale = cScale;
var_out.dScale = dScale;
var_out.D      = D;
var_out.E      = E;

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

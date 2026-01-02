function [runHist, sigma] = solver_socp_sGSinPALM(var, opts, model)
%% A sGS-based inPALM for solving the SOCP reformulation of Dynamic Optimal Transport:
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
sGSits      = 1;
updateRule  = [
    1.5, 1.20;
      2, 1.26;
    2.5, 1.28;
   3.33, 1.32;
      5, 1.35;
     10, 1.40;
];

% Scaling
cScale   = var.cScale;
dScale   = var.dScale;
D        = var.D;
E        = var.E;
scaleBF  = E / D;
scaleD   = E / dScale;
scaleLap = D^2;
use_feasOrg = false;
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

% Adjust sigma for sGS iteration
hist              = 19;
victory           = 12;
initialSigmaScale = 1.10;
stablePhase       = false;

%% initialization
% discrete model
nx      = model.nx;
ny      = model.ny;
nt      = model.nt;
h       = 1 / (nx*ny*nt);
A       = model.grad; % Grad
c       = model.c;

% iterative variable (detach to avoid struct copies; var is a handle)
phi     = var.phi;   var.phi   = [];
q       = var.q;     var.q     = [];
z       = var.z;     var.z     = [];
alpha   = var.alpha; var.alpha = [];
beta    = var.beta;  var.beta  = [];

% preprocessing
diagQInv = 1 ./ oper_q(ny, nx, nt, D, E);

%% Iteration
norm_c = model.normc;
norm_d = model.normd;
alpha  = alpha / sigma;
beta   = beta  / sigma;
c      = c     / sigma;
sigmaScale = 1;

sigma_adjust_it_gap = max(1, power(nt*nx*ny, 1/3) / 33);
sigma_adjust_val_gap = 0.95;
sgs_superior_yes = false;
tol_sgs_blocks = 5 * tol;

% Running history
kktConst = 1;
runHist  = {};
runHistItems  = 0;
runHist.kkt   = Inf(maxit, 7);
runHist.time  = Inf(maxit, 1);
runHist.iter  = Inf(maxit, 1);
runHist.pdGap = Inf(maxit, 1);
FeasRatio     = Inf(maxit, 1);
    
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
mexBFd(z2, q, nt, nx, ny, scaleBF, scaleD);

% Initial var
phi = phi - integralL2(phi, h);

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

        %   scale param
        sigma   = sigma  * (cScale2 / dScale2);
        c       = c      * dScale2 / cScale2^2;
        norm_c  = norm_c / cScale2;
        norm_d  = norm_d / dScale2;

        %   scale var
        alpha   = alpha  * dScale2 / cScale2^2;
        beta    = beta   * dScale2 / cScale2^2;
        phi     = phi    / dScale2;
        q       = q      / dScale2;

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
    mexsGS(phi, A' * (q - alpha) + c, 0, scaleLap, nt, nx, ny, sGSits);
    time_lineq  = time_lineq + toc(clock_lineq);
    
    clock_kkt = tic();
    adjustSigmaYes = IfAdjustSigma(it, lastSigmaIt, sigma_adjust_it_gap);
    check_kkt_yes = checkSByS || adjustSigmaYes || (it == maxit) || (toc(clock_total) > time_limit);
    if check_kkt_yes
        % Error of sGS blocks
        tmp_resi_sGS = A' * (A * phi - q + alpha) - c;
        % resi_sGS_1 = normL2(tmp_resi_sGS(1:2:end), h);
        % resi_sGS_2 = normL2(tmp_resi_sGS(2:2:end), h);
        resi_sGS_blocks = normL2(tmp_resi_sGS(1:2:end), h);
    end
    time_kkt = time_kkt + toc(clock_kkt);
    
    % step z
    clock_proj = tic();
    mexProjSoc(z, z2 - beta);
    time_proj = time_proj + toc(clock_proj);
    
    % step q
    clock_q = tic();
    tmp_q = A * phi;
    mexBFdConj(q2, z + beta, nt, nx, ny, scaleBF);
    q = ( tmp_q + alpha + q2 ) .* diagQInv;
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
    if check_kkt_yes
        % Precomputation
        %   temp
        mexBFdConj(q2, beta, nt, nx, ny, scaleBF);
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
        mexBFd(z2, q, nt, nx, ny, scaleBF, scaleD);
        
        [dotcomplem, normRho, norm_rhoFq, mRhoB, normM, normRhoB] = compute_kkt_dot_complement(q, alpha, z2, sigma, h, nt, nx, ny, var.qInd, cScale, dScale, D, E);
        
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
            primFea1  / (kktConst + norm_Aphi + norm_q), ...
            primFea2  / (kktConst + norm_d), ...
            dualFea1  / (kktConst + norm_c), ...
            complem   / (kktConst + norm_z + norm_beta), ...
            dualFea2  / (kktConst + norm_FBbeta + norm_alpha)
        ];
        FeasRatio(it) = max(KKTResi([1,2])) / max(KKTResi([3,5]));

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
            fprintf("Sigma: %.2E. PrimVal: %.2E. DualVal: %.2E. Prim-Dual gap: %.2E. PrimFeas: %.2E. DualFeas: %.2E. Complement: %.2E.\n", ...
                sigma, priVal, dualVal, pdGap, max(KKTResiOrg([1,2])), max(KKTResiOrg([3,5])), KKTResiOrg(4));
        end

        if (printDotYes)
            fprintf("PrimVal: %.2E. DualVal: %.2E. Prim-Dual gap: %.2E. DualFeas: %.2E. PrimFeas: %.2E. Compleme: %.2E. PrimDualFeas: %.2E.\n", ...
                priVal, dualVal, pdGap, KKTResiOrg(1), KKTResiOrg(3), KKTResiOrg(6), KKTResiOrg(7));
        end

        % Stop criterion
        error = max(KKTResiOrg(stopCondition));
        if error < tol || toc(clock_total) > time_limit
            break;
        end

        % Use original KKT
        if ~use_feasOrg && max(KKTResi) < tol_feasOrg
            use_feasOrg = true;
        end
    
        % Update Lagrangian parameter
        kkt_sgs_blocks = sqrt(normL2(A' * resi_alpha, h)^2 + (dualFea1 / sigma)^2);
        sgs_superior_yes = resi_sGS_blocks < sigma_adjust_val_gap * kkt_sgs_blocks;

        if (printYes)
            fprintf("(Residual, KKT error) of sGS blocks: (%.2E, %.2E) with superior ratio: %.2E.\n", resi_sGS_blocks, kkt_sgs_blocks, resi_sGS_blocks / kkt_sgs_blocks);
        end

        if adjustSigmaYes
            lastSigmaIt = it;

            feasRatioHist = FeasRatio(max(1,it-hist) : it);
            meanFeasRatio = mean(feasRatioHist);
            primWinTimes  = sum(feasRatioHist < 1);
            dualWinTimes  = sum(feasRatioHist > 1);

            adjust_sigma_yes_2 = sgs_superior_yes || (error < tol_sgs_blocks) || ...
                (dualWinTimes >= victory) && (meanFeasRatio > 1);
            
            if adjust_sigma_yes_2
                if (it > 2500)
                    stablePhase = true;
                end
    
                if (primWinTimes >= victory) && (meanFeasRatio < 1) ...
                        || (dualWinTimes >= victory) && (meanFeasRatio > 1)
                    if stablePhase
                        [sigma, factor] = adjust_lagrangianParam(sigma, meanFeasRatio, updateRule);
                    else
                        if (meanFeasRatio < 1)
                            factor = 1 / initialSigmaScale ;
                        elseif (meanFeasRatio > 1)
                            factor = initialSigmaScale;
                        end
    
                        sigma  = sigma * factor;
                    end
    
                    if factor ~= 1
                        alpha  = alpha / factor;
                        beta   = beta  / factor;
                        c      = c     / factor;
                    end
                end
            end
        end

        % Information for rescaling
        if (rescale > 0)
            maxFeas = max(KKTResi);
            relGap  = pdGap;
        end
    elseif sgs_superior_yes
        % tmp_q     = A * phi;
        % norm_q    = normL2(q, h);
        % norm_Aphi = normL2(tmp_q, h);
        primFea1  = normL2(resi_alpha, h);
        dualFea1  = sigma * normL2(A'*alpha - c, h);

        if use_feasOrg
            relaPrimFeaDec = primFea1 / ( (kktConst * D / dScale + norm_Aphi + norm_q) * KKTResi(1));
            KKTResi([1,2]) = KKTResi([1,2]) * relaPrimFeaDec;
            KKTResi(3) = dualFea1 / (kktConst / cScale + norm_c);
        else
            relaPrimFeaDec = primFea1 / ( (kktConst + norm_Aphi + norm_q) * KKTResi(1));
            KKTResi([1,2]) = KKTResi([1,2]) * relaPrimFeaDec;
            KKTResi(3) = dualFea1 / (kktConst + norm_c);
        end

        FeasRatio(it) = max(KKTResi([1,2])) / max(KKTResi([3,5]));
    else
        FeasRatio(it) = FeasRatio(it-1);
    end
    time_kkt = time_kkt + toc(clock_kkt);
end
time_total = toc(clock_total);

%% output
var.name = 'Symmetric Gauss-seidel based inPALM';

% Iterative var
var.phi = phi;
var.q = q;
var.z = z;
var.alpha = sigma * alpha;
var.beta = sigma * beta;

% Time
times = [time_lineq, time_proj, time_q, time_multiplier, time_kkt, time_total, it];
names = {'Step_1_1_sGS', 'Step_1_2_ProjSOC', 'Step_2_Q_Step', 'Step_3_Multiplier', 'KKT', 'Total_Time', 'Iters'};
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

function flag = IfAdjustSigma(iter, last_iter, scale)
    %% Whether to adjust Lagrangian param
    passedIters = iter - last_iter;
    flag = 0;

    if nargin <= 2
        scale = 1;
    end

    iter = iter / scale;
    passedIters = passedIters / scale;
    
    if iter < 20 && passedIters >= 5
        flag = 1;
    elseif iter < 50 && passedIters >= 10
        flag = 1;
    elseif iter < 100 && passedIters >= 20
        flag = 1;
    elseif iter < 200 && passedIters >= 35
        flag = 1;
    elseif iter < 500 && passedIters >= 50
        flag = 1;
    elseif passedIters >= 100
        flag = 1;
    end
end

function [var_out, runHist, sigma] = solver_socp_accsGSADMM(var, opts, model)
%% An accelerated sGS-based ADMM for solving the SOCP reformulation of Dynamic Optimal Transport:
%       min <c, \phi> + \delta_{Q}(z)
%       s.t.    A \phi - q = 0,
%               z - B F q  = d.
% *************************************************************************
% Copyright (c) 2024 by
% Liang Chen, Youyicun Lin, and Yuxuan Zhou
% *************************************************************************

%% Params
if ~isfield(opts, 'restart')
    restart = 100;
else
    restart = opts.restart;
end

if ~isfield(opts, 'rho')
    stepRho = 2;
else
    stepRho = opts.rho;
end

if ~isfield(opts, 'theta')
    stepAlpha = 2;
else
    stepAlpha  = opts.theta;
end

if (stepAlpha == 2)
    HalpernYes = true; % Halpern iteration
else
    HalpernYes = false;
end

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
     20, 1.60;
     40, 1.80;
     50, 2.00;
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
    checkRescaleIters = 200;
    ratioThreshold    = 1.2;
    maxFeas           = Inf;
    relGap            = Inf;
end

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
At      = transpose(A); % -Div
c       = model.c;

sigma_adjust_it_gap = max(1, power(nt*nx*ny, 1/3) / 33);
sigma_adjust_val_gap = 0.95;
sgs_superior_yes = false;
tol_sgs_blocks = 5 * tol;

% iterative variable
phi     = var.phi;
q       = var.q;
z       = var.z;
alpha   = var.alpha;
beta    = var.beta;

% preprocessing
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
time_interp     = 0;

% Initial var
phi = phi - integralL2(phi, h);

% Preallocation
z2 = zeros((nt-1)*nx*ny, 10);
q2 = zeros((nt-1)*nx*ny + nt*(nx-1)*ny + nt*nx*(ny-1), 1);
projZBta = zeros((nt-1)*nx*ny, 10);
[phiOld, zOld, qOld, alphaOld, betaOld] = CopyVar(phi, z, q, alpha, beta);
k = 0;

% set initial anchor for Halpern
if HalpernYes
    [phi0, z0, q0, alpha0, beta0] = CopyVar(phi, z, q, alpha, beta);
end

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
        z       = z      / dScale2;

        %   record scaling factor
        dScale  = dScale2 * dScale;
        cScale  = cScale2 * cScale;
        scaleD  = E / dScale;
        sigmaScale = sigmaScale * (cScale2 / dScale2);

        % restart
        k = 0;
        [phiOld, zOld, qOld, alphaOld, betaOld] = CopyVar(phi, z, q, alpha, beta);
        if HalpernYes % reset anchor
            [phi0, z0, q0, alpha0, beta0] = CopyVar(phi, z, q, alpha, beta);
        end

        rescale = rescale + 1;
    end

    % step q
    clock_q = tic();
    mexBFdConj(q2, z + beta, nt, nx, ny, scaleBF);
    tmp_q = A * phi;
    q = (tmp_q + alpha + q2) .* diagQInv;
    time_q = time_q + toc(clock_q);

    % step alpha, beta
    clock_multiplier = tic();
    mexBFd(z2, q, nt, nx, ny, scaleBF, scaleD);
    alpha = alpha + tmp_q - q;
    beta  = beta  + z - z2;
    time_multiplier = time_multiplier + toc(clock_multiplier);

    % step phi
    clock_lineq = tic();
    mexsGS(phi, At * (q - alpha) + c, 0, scaleLap, nt, nx, ny, sGSits);
    time_lineq  = time_lineq + toc(clock_lineq);

    clock_kkt = tic();
    adjustSigmaYes = IfAdjustSigma(it, lastSigmaIt, sigma_adjust_it_gap);
    check_kkt_yes = checkSByS || adjustSigmaYes || (it == maxit) || (toc(clock_total) > time_limit);
    if check_kkt_yes
        % Error of sGS blocks
        tmp_resi_sGS = At * (A * phi - q + alpha) - c;
        % resi_sGS_1 = normL2(tmp_resi_sGS(1:2:end), h);
        % resi_sGS_2 = normL2(tmp_resi_sGS(2:2:end), h);
        resi_sGS_blocks = normL2(tmp_resi_sGS(1:2:end), h);
    end
    time_kkt = time_kkt + toc(clock_kkt);

    % step z
    clock_proj = tic();
    mexProjSoc(z, z2 - beta);
    time_proj = time_proj + toc(clock_proj);
    
    % step interpolation
    clock_interp = tic();
    if HalpernYes
        % Halpern Iteration
        c1    =   1   / (k+2);
        c2    = (k+1) / (k+2);
        phi   = c1 * phi0   + c2 * ( (1-stepRho) * phiOld   + stepRho * phi  );
        z     = c1 * z0     + c2 * ( (1-stepRho) * zOld     + stepRho * z    );
        q     = c1 * q0     + c2 * ( (1-stepRho) * qOld     + stepRho * q    );
        alpha = c1 * alpha0 + c2 * ( (1-stepRho) * alphaOld + stepRho * alpha);
        beta  = c1 * beta0  + c2 * ( (1-stepRho) * betaOld  + stepRho * beta );

        k = k + 1;
        [phiOld, zOld, qOld, alphaOld, betaOld] = CopyVar(phi, z, q, alpha, beta);
    
        % restart
        if k >= restart
            k = 0;
            [phi0, z0, q0, alpha0, beta0] = CopyVar(phi, z, q, alpha, beta);
        end
    else
        phiHat   = (1-stepRho) * phiOld   + stepRho * phi;
        zHat     = (1-stepRho) * zOld     + stepRho * z;
        qHat     = (1-stepRho) * qOld     + stepRho * q;
        alphaHat = (1-stepRho) * alphaOld + stepRho * alpha;
        betaHat  = (1-stepRho) * betaOld  + stepRho * beta;
    
        if (k == 0)
            c1 = stepAlpha / (2 * (k + stepAlpha));
            phi   = (1 - c1) * phiOld   + c1 * phiHat  ;
            z     = (1 - c1) * zOld     + c1 * zHat    ;
            q     = (1 - c1) * qOld     + c1 * qHat    ;
            alpha = (1 - c1) * alphaOld + c1 * alphaHat;
            beta  = (1 - c1) * betaOld  + c1 * betaHat ;
        else
            c1    = stepAlpha / (2 * (k + stepAlpha));
            c2    = k / (k + stepAlpha);
            phi   = (1 - c1) * phiOld   + (c1 + c2) * phiHat   - c2 * phiHatOld  ;
            z     = (1 - c1) * zOld     + (c1 + c2) * zHat     - c2 * zHatOld    ;
            q     = (1 - c1) * qOld     + (c1 + c2) * qHat     - c2 * qHatOld    ;
            alpha = (1 - c1) * alphaOld + (c1 + c2) * alphaHat - c2 * alphaHatOld;
            beta  = (1 - c1) * betaOld  + (c1 + c2) * betaHat  - c2 * betaHatOld ;
        end

        k = k + 1;
        [phiOld, zOld, qOld, alphaOld, betaOld] = CopyVar(phi, z, q, alpha, beta);

        % restart
        if k >= restart
            k = 0;
        else
            [phiHatOld, zHatOld, qHatOld, alphaHatOld, betaHatOld] = deal(phiHat, zHat, qHat, alphaHat, betaHat);
        end
    end
    time_interp = time_interp + toc(clock_interp);
    
    % kkt
    clock_kkt = tic();
    if check_kkt_yes
        % Precomputation
        %   temp
        mexBFdConj(q2, beta, nt, nx, ny, scaleBF);
        mexBFd(z2, q, nt, nx, ny, scaleBF, scaleD);
        mexProjSoc(projZBta, z - sigma * beta);
        tmp_q = A * phi;
        resi_alpha = tmp_q - q;
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
        primFea2   = FnormL2(z - z2, h);
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

        % Print iterative message
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
        kkt_sgs_blocks = sqrt(normL2(At * resi_alpha, h)^2 + (dualFea1 / sigma)^2);
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
                if (it > 1500)
                    stablePhase = true;
                end
    
                % Adjust sigma
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
                        [alpha, alphaOld] = deal(alpha / factor);
                        [beta,  betaOld ] = deal(beta / factor);
                        c                 = c / factor;
    
                        % restart
                        k = 0;
                        if HalpernYes
                            [phi0, z0, q0, alpha0, beta0] = CopyVar(phi, z, q, alpha, beta);
                        end
                    end
                end
            end
        end

        % Infomation for rescaling
        if (rescale > 0)
            maxFeas = max(KKTResi);
            relGap  = pdGap;
        end
    elseif sgs_superior_yes
        % Update a part of primFeas and dualFeas
        % norm_q    = normL2(q, h);
        % norm_Aphi = normL2(tmp_q, h);
        primFea1  = normL2(A * phi - q, h);
        dualFea1  = sigma * normL2(At*alpha - c, h);

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
var_out = var;
var_out.name = 'Accelerated symmetric Gauss-Seidel based ADMM';

% Iterative var
var_out.phi = phi;
var_out.q = q;
var_out.z = z;
var_out.alpha = sigma * alpha;
var_out.beta = sigma * beta;

% Time
times = [time_lineq, time_proj, time_multiplier, time_q, time_interp, time_kkt, time_total, it];
names = {'Step_1_1_sGS', 'Step_1_2_ProjSOC', 'Step_2_Multiplier' , 'Step_3_Q_Step', 'Step_4_Interp', 'KKT', 'Total_Time', 'Iters'};
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

function [phiCopy, zCopy, qCopy, alphaCopy, betaCopy] = CopyVar(phi, z, q, alpha, beta)
    %% Copy variable
    [qCopy, alphaCopy, betaCopy] = deal(q, alpha, beta);
    phiCopy = phi(1:end);
    zCopy   = reshape(z(1:end), size(z));
end
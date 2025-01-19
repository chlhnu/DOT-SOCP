function [varR, modelR] = jump_nextLevel(var, model, rho0, rho1, nt)
%% Jump into next level in the multilevel strategy

% Interpolate variables to finer grid
varR = interpolate(var, model);

% Initialize refined discrete model
[var_init, modelR] = initialize(rho0, rho1, nt);
varR.qInd = var_init.qInd;
varR.z    = var_init.z;

% Preprocess variables
varR.q = modelR.grad * varR.phi;
varR.alpha = var_init.alpha; % preallocation
mexBFdConj1d(varR.alpha, - varR.beta, modelR.nt, modelR.nx, 1);

end


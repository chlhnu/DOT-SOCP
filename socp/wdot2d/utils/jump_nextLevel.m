function [varR, modelR] = jump_nextLevel(var, model, rho0, rho1, weight, nt)
%% Jump to the next level in the multilevel scheme

varR = interpolate(var, model);

[var_init, modelR] = initialize(rho0, rho1, nt);
varR.qInd = var_init.qInd;
varR.z    = var_init.z;
modelR.weight = weight;

varR.q = (modelR.grad * varR.phi) ./ weight;
varR.alpha = var_init.alpha; % preallocation
mexBFdConj(varR.alpha, - varR.beta, modelR.nt, modelR.nx, modelR.ny, 1);
varR.alpha = varR.alpha ./ weight;

end


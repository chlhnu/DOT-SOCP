function [rho, Ex] = recover_RhoE(var, model)
%% Recover (rho, Ex) from alpha

nt    = model.nt;
nx    = model.nx;
rho0  = model.rho0;
rho1  = model.rho1;

qInd  = var.qInd;
alpha = var.alpha;

% rho
rho = reshape(alpha(1:qInd.bx-1), nx, nt-1);
rho = cat(2, rho0, (rho(:,1:end-1) + rho(:,2:end))/2, rho1);

% Ex
Ex = reshape(alpha(qInd.bx : end), nx-1, nt);
Ex(:,[1,end]) = 2 * Ex(:,[1,end]);
Ex = cat(1, zeros(1,nt), (Ex(1:end-1,:) + Ex(2:end,:))/2, zeros(1,nt));

end

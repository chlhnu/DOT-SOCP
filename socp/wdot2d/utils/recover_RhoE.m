function [rho, Ex, Ey] = recover_RhoE(var, model)
%% Recover (rho, Ex, Ey) from alpha

nt    = model.nt;
nx    = model.nx;
ny    = model.ny;
rho0  = model.rho0;
rho1  = model.rho1;

qInd  = var.qInd;
alpha = model.weight .* var.alpha;

% rho
rho = reshape(alpha(1:qInd.bx-1), ny, nx, nt-1);
rho = cat(3, rho0, (rho(:,:,1:end-1) + rho(:,:,2:end))/2, rho1);

% Ex
Ex = reshape(alpha(qInd.bx : qInd.by-1), ny, nx-1, nt);
Ex(:,:,[1,end]) = 2 * Ex(:,:,[1,end]);
Ex = cat(2, zeros(ny,1,nt), (Ex(:,1:end-1,:) + Ex(:,2:end,:))/2, zeros(ny,1,nt));

% Ey
Ey = reshape(alpha(qInd.by : end), ny-1, nx, nt);
Ey(:,:,[1,end]) = 2 * Ey(:,:,[1,end]);
Ey = cat(1, zeros(1,nx,nt), (Ey(1:end-1,:,:) + Ey(2:end,:,:))/2, zeros(1,nx,nt));

end

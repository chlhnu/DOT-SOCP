function [q0, bx, by] = recover_q(var, model)
%% Recover (q0, bx, by) from staggered grid to time-staggered grid

nt    = model.nt;
nx    = model.nx;
ny    = model.ny;

qInd  = var.qInd;
q = var.q;

% Recover q0 component
q0 = reshape(q(1:qInd.bx-1), ny, nx, nt-1);

% Recover bx component with boundary conditions
bx = reshape(q(qInd.bx : qInd.by-1), ny, nx-1, nt);
bx = cat(2, zeros(ny,1,nt), (bx(:,1:end-1,:) + bx(:,2:end,:))/2, zeros(ny,1,nt));
bx = (bx(:, :, 1:end-1) + bx(:, :, 2:end)) / 2;

% Recover by component with boundary conditions
by = reshape(q(qInd.by : end), ny-1, nx, nt);
by = cat(1, zeros(1,nx,nt), (by(1:end-1,:,:) + by(2:end,:,:))/2, zeros(1,nx,nt));
by = (by(:, :, 1:end-1) + by(:, :, 2:end)) / 2;

end

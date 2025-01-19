function [q0, bx] = recover_q(var, model)
%% Recover (q0, bx)

nt    = model.nt;
nx    = model.nx;
qInd  = var.qInd;
q     = var.q;

% Recover q0 component
q0 = reshape(q(1:qInd.bx-1), nx, nt-1);

% Recover bx component with boundary conditions
bx = reshape(q(qInd.bx : end), nx-1, nt);
bx = cat(1, zeros(1,nt), (bx(1:end-1,:) + bx(2:end,:))/2, zeros(1,nt));
bx = (bx(:, 1:end-1) + bx(:, 2:end)) / 2;

end

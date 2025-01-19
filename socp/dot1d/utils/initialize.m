function [Var, Model] = initialize(rho0, rho1, nt)
%% Initialization for DOT-SOCP (1 dimension)
Var   = {};
Model = {};

% params
nx = length(rho0);
n = nx * nt;

% Rho
Model.rho0 = rho0;
Model.rho1 = rho1;

rho0 = rho0(:);
rho1 = rho1(:);

%% params
qInd = {};
qInd.bx  = (nt-1) * nx + 1;

% stepsize
ht = 1/(nt-1);
hx = 1/(nx-1);

% output
Model.nx = nx;
Model.nt = nt;
Var.qInd = qInd;

%% Poisson
% Grad
gradt = gene_Dt(nt, nx, ht);
gradx = gene_Dx(nt, nx, hx);
grad = [gradt; gradx];

% -c
c = zeros(n, 1);
c(1:nx) = - rho0 / ht;
c(end-nx+1:end) = rho1 / ht;

% output
Model.grad = grad; % grad
Model.c    = c;    % -c

%% Initialize var
% phi
xx = 0 : hx : 1;
Var.phi      = 1/2 * xx.^2;
Var.phi      = reshape(repmat(Var.phi(:), [1, nt]), [], 1);

% z, beta
lenA         = (nt-1) * nx;
Var.z        = zeros(lenA, 6);
Var.beta     = zeros(lenA, 6);

% q, alpha
Var.q        = zeros(size(grad, 1), 1);
Var.alpha    = zeros(size(grad, 1), 1);

end

function gradt = gene_Dt(nt, nx, ht)
    tmp = repmat(1/ht, [nt, 1]);
    Dt = spdiags([-tmp, tmp], [0, 1], nt-1, nt);
    Ix = speye(nx);
    gradt = kron(Dt, Ix);
end

function gradx = gene_Dx(nt, nx, hx)
    It = speye(nt);
    tmp = repmat(1/hx, [nx, 1]);
    Dx = spdiags([-tmp, tmp], [0, 1], nx-1, nx);
    gradx = kron(It, Dx);
end

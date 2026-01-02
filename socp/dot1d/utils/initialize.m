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
Model.grad = [
    gene_Dt(nt, nx, ht);
    gene_Dx(nt, nx, hx)
];

% -c
Model.c = zeros(n, 1);
Model.c(1:nx) = - rho0(:) / ht;
Model.c(end-nx+1:end) = rho1(:) / ht;

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
Var.q        = zeros(size(Model.grad, 1), 1);
Var.alpha    = zeros(size(Model.grad, 1), 1);

%% To handle
Var = VarHandle(Var);
Model = ModelHandle(Model);

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

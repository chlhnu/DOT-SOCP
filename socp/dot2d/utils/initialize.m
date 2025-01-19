function [var, model] = initialize(rho0, rho1, nt)
%% Initialization for DOT-SOCP

var   = {};
model = {};

% params
nx = size(rho0, 2);
ny = size(rho0, 1);
n = nx*ny*nt;

% Rho
model.rho0 = rho0;
model.rho1 = rho1;

% vectorization
rho0 = rho0(:);
rho1 = rho1(:);

%% params
% index of (bx, by)
qInd = {};
qInd.bx  = (nt-1) * nx * ny + 1;
qInd.by  = nt * (nx-1) * ny + qInd.bx;

% stepsize
ht = 1/(nt-1);
hx = 1/(nx-1);
hy = 1/(ny-1);

% output
model.nx = nx;
model.ny = ny;
model.nt = nt;
var.qInd = qInd;

%% Poisson
% Grad
gradt = gene_Dt(nt, nx, ny, ht);
gradx = gene_Dx(nt, nx, ny, hx);
grady = gene_Dy(nt, nx, ny, hy);
grad = [gradt; gradx; grady];

% -c
c = zeros(n, 1);
c(1:nx*ny) = -rho0 / ht;
c(end-nx*ny+1:end) = rho1 / ht;

% output
model.grad = grad; % grad
model.c    = c;    % -c

%% Initlization
% phi
[xx, yy]     = meshgrid(0 : hx : 1, 0 : hy : 1);
var.phi      = 1/2 * (xx.^2 + yy.^2);
var.phi      = reshape(repmat(var.phi, [1, 1, nt]), [], 1);

% z, beta
lenA         = (nt-1) * nx * ny;
var.z        = zeros(lenA, 10);
var.beta     = zeros(lenA, 10);

% q, alpha
var.q        = zeros(size(grad, 1), 1);
var.alpha    = zeros(size(grad, 1), 1);

end

function gradt = gene_Dt(nt, nx, ny, ht)
    tmp = repmat(1/ht, [nt, 1]);
    Dt = spdiags([-tmp, tmp], [0, 1], nt-1, nt);
    Ixy = speye(nx*ny);
    gradt = kron(Dt, Ixy);
end

function gradx = gene_Dx(nt, nx, ny, hx)
    It = speye(nt);
    Iy = speye(ny);
    tmp = repmat(1/hx, [nx, 1]);
    Dx = spdiags([-tmp, tmp], [0, 1], nx-1, nx);
    gradx = kron(kron(It, Dx), Iy);
end

function grady = gene_Dy(nt, nx, ny, hy)
    Itx = speye(nt*nx);
    tmp = repmat(1/hy, [ny, 1]);
    Dy = spdiags([-tmp, tmp], [0, 1], ny-1, ny);
    grady = kron(Itx, Dy);
end

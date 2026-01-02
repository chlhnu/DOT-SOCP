function [rho0, rho1] = gene_exampleCircle2(nx, ny)
%% Generate circular densities (used in barrier of rectangle)

hx = 1/(nx-1);
hy = 1/(ny-1);
[xx, yy] = meshgrid(0:hx:1, 0:hy:1);

scale = 40;
r1 = 5 / scale;
r2 = 4 / scale;
r3 = 3 / scale;

% rho0
rho0 = zeros(ny, nx);

center0 = [r1/2 + 0.1, 0.475];
mask0 = (xx-center0(1)).^2 + (yy-center0(2)).^2 < r1^2;
rho0(mask0) = 1;

% rho1
rho1 = zeros(ny, nx);

center1 = [r2/2 + 0.1, 0.95 - r2];
mask11 = (xx-center1(1)).^2 + (yy-center1(2)).^2 < r2^2;
rho1(mask11) = 1;

center2 = [r3/2 + 0.1, r3 + 0.05];
mask12 = (xx-center2(1)).^2 + (yy-center2(2)).^2 < r3^2;
rho1(mask12) = 1;

end


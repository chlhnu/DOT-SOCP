function [rho0, rho1] = gene_exampleCircle(nx, ny)
%% Generate circular densities

% grid, x/y is the 2/1-th index
hx = 1/(nx-1);
hy = 1/(ny-1);
[xx, yy] = meshgrid(0:hx:1, 0:hy:1);

% radius and center of circle
r1 = 0.25;
r2 = 0.25;
center1 = [0.25, 0.75];
center2 = [0.75, 0.25];

% rho
rho0 = zeros(ny, nx);
mask0 = (xx-center1(1)).^2 + (yy-center1(2)).^2 < r1^2;
rho0(mask0) = 1;

rho1 = zeros(ny, nx);
mask1 = (xx-center2(1)).^2 + (yy-center2(2)).^2 < r2^2;
rho1(mask1) = 1;

end


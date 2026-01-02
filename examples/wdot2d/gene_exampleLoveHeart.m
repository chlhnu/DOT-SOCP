function [rho0, rho1] = gene_exampleLoveHeart(nx, ny)
%% Generate Gaussian densities for Example of Love Heart

% center, radius
center1 = [0.7, 0.3];
center2 = [0.345, 0.625];
r1 = 0.09;
r2 = 0.09;

% sigma1, sigma2 in (0, 0.2)
sigma1  = r1 / 3;
sigma2  = r2 / 3;

% Gaussian
[Y, X] = meshgrid(linspace(0, 1, nx), linspace(0, 1, ny));
gaussian = @(a, b, sigma) exp( - ( (X-a).^2 + (Y-b).^2 ) / (2*sigma^2) );

% rho0
rho0 = gaussian(center1(1), center1(2), sigma1);
mask = (X-center1(1)).^2 + (Y-center1(2)).^2 > r1^2;
rho0(mask) = 0;

rho1 = gaussian(center2(1), center2(2), sigma2);
mask = (X-center2(1)).^2 + (Y-center2(2)).^2 > r2^2;
rho1(mask) = 0;

end


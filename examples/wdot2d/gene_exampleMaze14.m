function [rho0, rho1] = gene_exampleMaze14(nx, ny)
%% Generate Gaussian densities, similar to an example from [Optimal Transport with Proximal Splitting. SIAM Journal on Imaging Sciences, 2014.]

% center, radius
center1 = [0.075, 0.075];
center2 = [0.925, 0.925];
r1 = 0.075;
r2 = 0.075;

% sigma1, sigma2 in (0, 0.2)
sigma1  = r1 / 2;
sigma2  = r2 / 2;

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


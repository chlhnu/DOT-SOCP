function [rho0, rho1] = gene_example4(nx, ny)
%% Generate densities for Example 5.4

% grid
[Y, X] = meshgrid(linspace(0, 1, nx), linspace(0, 1, ny));

% rho0
center = [0.5, 0.5];
rho0 = (X-center(1)).^4 + (Y-center(2)).^4;

% rho1
mu1    = 0.25;
mu2    = 1 - mu1;
sigma2 = 0.05;
gaussian = @(a, b, sigma) exp( - ( (X-a).^2 + (Y-b).^2 ) / (2*sigma^2) );
rho1 = gaussian(mu1, mu1, sigma2) + gaussian(mu1, mu2, sigma2) + gaussian(mu2, mu1, sigma2) + gaussian(mu2, mu2, sigma2);

end


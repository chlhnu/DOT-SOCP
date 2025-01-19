function [rho0, rho1] = gene_example3(nx, ny)
%% Example 5.3

[a1, a2] = deal(3, 5);

% mu1, mu2 in (0,1) and mu1 + mu2 == 1
mu1 = 0.25;
mu2 = 1 - mu1;

% sigma in (0, 0.2)
sigma  = 0.05;

% grid
[Y, X] = meshgrid(linspace(0, 1, nx), linspace(0, 1, ny));

% rho0
laplacian = @(a, b) exp( exp( - a1 * abs(X-a) - a2 * abs(Y-b) ) );
rho0 = laplacian(mu1, mu1);

% rho1
gaussian = @(a, b, sigma) exp( - ( (X-a).^2 + (Y-b).^2 ) / (2*sigma^2) );
rho1 = gaussian(mu1, mu1, sigma) + gaussian(mu1, mu2, sigma) + gaussian(mu2, mu1, sigma) + gaussian(mu2, mu2, sigma);

end


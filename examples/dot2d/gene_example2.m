function [rho0, rho1] = gene_example2(nx, ny)
%% Example 5.2

% mu1, mu2 in (0,1) and mu1 + mu2 == 1
mu1 = 0.25;
mu2 = 1 - mu1;

% sigma1, sigma2 in (0, 0.2)
sigma1  = 0.1;
sigma2  = 0.05; % sigma1 approx 4 * sigma2^2

% grid
[Y, X] = meshgrid(linspace(0, 1, nx), linspace(0, 1, ny));

% rho
gaussian = @(a, b, sigma) exp( - ( (X-a).^2 + (Y-b).^2 ) / (2*sigma^2) );
rho0 = gaussian(mu1, mu1, sigma1);
rho1 = gaussian(mu1, mu1, sigma2) + gaussian(mu1, mu2, sigma2) + gaussian(mu2, mu1, sigma2) + gaussian(mu2, mu2, sigma2);

end


function [rho0, rho1] = gene_example_gaussian(nx)
%% Generate Gaussian densities

% mu1, mu2 in (0, 1)
mu1 = 0.3;
mu2 = 0.7;

% sigma
sigma1 = 0.01;
sigma2 = sigma1 / 4;

% Gaussian
Normal = @(x, meanvalue, sigmainv) sqrt(sigmainv) / (2*pi) * ...
    exp( -0.5 * (sigmainv * (x - meanvalue).^2) );

% grid
x = reshape(linspace(0,1,nx), nx, 1);

% rho
rho0 = Normal(x, mu1, 1 / sigma1);
rho1 = Normal(x, mu2, 1 / sigma2);

end

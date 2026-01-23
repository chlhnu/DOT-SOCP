function [rho0, rho1] = gene_example1(nx, ny)
%% Example 5.1

% mu1, mu2 in (0, 1)
mu1 = 0.25;
mu2 = 0.75;

% sigma
sigma = 0.05;

% Gaussian
meanvalue0 = [mu1,mu2]; sigma0 = [sigma,0;0,sigma];
meanvalue1 = [mu2,mu1]; sigma1 = [sigma,0;0,sigma];
Normal = @(x,y,meanvalue,sigmainv) sqrt(det(sigmainv))/(2*pi) * ...
    exp(-0.5* (sigmainv(1,1)*(x-meanvalue(1)).^2 + ...
               sigmainv(1,2)*(x-meanvalue(1)).*(y-meanvalue(2)) + ...
               sigmainv(2,2)*(y-meanvalue(2)).^2));

% grid
x = repmat(reshape(linspace(0,1,nx), nx, 1), 1, ny);
y = repmat(reshape(linspace(0,1,ny), 1, ny), nx, 1);

% rho
rho0 = Normal(x,y,meanvalue0,inv(sigma0));
rho1 = Normal(x,y,meanvalue1,inv(sigma1));

end

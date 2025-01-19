function [rho0, rho1] = gene_example7(nx, ny)
%% Example 5.7

% Four dirac points
% diracX = [0.25, 0.25, 0.75, 0.75];
% diracY = [0.25, 0.75, 0.25, 0.75];

% A row and a column dirac points
% diracX = [0.1, 0.3, 0.5, 0.7, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9];
% diracY = [0.9, 0.9, 0.9, 0.9, 0.9, 0.1, 0.3, 0.5, 0.7, 0.9];

% Random dirac points
%   Generator
% N = 30;
% origin = [0.5, 0.5];
% r = .48 * rand(1, N);
% theta = (2*pi) * rand(1, N);
% diracX = origin(1) + r .* cos(theta);
% diracY = origin(2) + r .* sin(theta);
%   An instance from above generator
diracX = [0.8323    0.5339    0.4031    0.6536    0.8200    0.4918    0.5108    0.6082    0.4633    0.1500    0.7227    0.4967    0.5318    0.6625    0.4309    0.1076    0.3052    0.4113    0.4955    0.4485    0.5031    0.7529    0.4723    0.3668    0.4848    0.5474    0.3867    0.3192    0.0676    0.2382];
diracY = [0.4477    0.6033    0.4264    0.5378    0.8026    0.7535    0.3472    0.2628    0.4023    0.4676    0.4535    0.5105    0.5903    0.6705    0.5134    0.4471    0.6960    0.5068    0.5040    0.5468    0.2641    0.1783    0.2195    0.3484    0.5056    0.3925    0.4511    0.2659    0.4157    0.8016];

%% Main
% grid, x/y is the 2/1-th index
hx = 1/(nx-1);
hy = 1/(ny-1);
[xx, yy] = meshgrid(0:hx:1, 0:hy:1);

% Rho0 - Gauss
[mux, muy] = deal(0.5, 0.5);
sigma = 0.1;
gaussian = @(a, b, sigma) exp( - ( (xx-a).^2 + (yy-b).^2 ) / (2*sigma^2) );
rho0 = gaussian(mux, muy, sigma);

% Rho1 - Dirac
rho1 = zeros(ny, nx);
diracXIndex = max(1, min(nx, round(diracX / hx)));
diracYIndex = max(1, min(nx, round(diracY / hy)));
for i = 1 : length(diracX)
    rho1(diracXIndex(i), diracYIndex(i)) = 1;
end

end


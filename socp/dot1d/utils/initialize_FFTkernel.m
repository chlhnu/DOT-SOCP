function kernel = initialize_FFTkernel(nt, nx, epsilon)
%% Get 2D discrete Fourier cosine transform kernel for Poisson/Helmholtz equation
% Optional input
%   epsilon, Helmholtz equation - u_tt - u_xx + epsilon * u = f, where epsilon >= 0

CT = (2 * (nt-1)^2) * (1 - cos(pi * (0:1:nt-1) / nt));
CX = (2 * (nx-1)^2) * (1 - cos(pi * (0:1:nx-1) / nx));

if nargin == 2 || epsilon == 0
    kernel = repmat(reshape(CX, nx, 1), [1, nt]) + ...
             repmat(reshape(CT, 1, nt), [nx, 1]);
    
    kernel(kernel(:) == 0) = 1;
elseif nargin == 3
    kernel = repmat(reshape(CX, nx, 1), [1, nt]) + ...
             repmat(reshape(CT, 1, nt), [nx, 1]) + ...
             epsilon;
else
    error("initialize_FFTkernel.m: Invalid input");
end

end


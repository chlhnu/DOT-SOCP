function kernel = initialize_FFTkernel(nt, nx, ny, epsilon)
%% Get 3D discrete Fourier cosine transform kernel for Poisson/Helmholtz equation
% Optional input:
%   epsilon - Parameter for Helmholtz equation: u_tt - u_xx - u_yy + epsilon*u = f, where epsilon >= 0

CT = (2 * (nt-1)^2) * (1 - cos(pi * (0:1:nt-1) / nt));
CX = (2 * (nx-1)^2) * (1 - cos(pi * (0:1:nx-1) / nx));
CY = (2 * (ny-1)^2) * (1 - cos(pi * (0:1:ny-1) / ny));

if nargin == 3
    kernel = repmat(reshape(CY, ny, 1, 1), [1, nx, nt]) + ...
             repmat(reshape(CX, 1, nx, 1), [ny, 1, nt]) + ...
             repmat(reshape(CT, 1, 1, nt), [ny, nx, 1]);
         
    kernel(kernel(:) == 0) = 1;
elseif nargin == 4
    kernel = repmat(reshape(CY, ny, 1, 1), [1, nx, nt]) + ...
             repmat(reshape(CX, 1, nx, 1), [ny, 1, nt]) + ...
             repmat(reshape(CT, 1, 1, nt), [ny, nx, 1]) + ...
             epsilon;
end

end

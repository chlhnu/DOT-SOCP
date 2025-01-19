function flag = check_massConservation(rho, tol, varargin)
%% Check mass conservation error in each time layer

p = inputParser;
addParameter(p, 'Show', false);
parse(p, varargin{:});
Show = p.Results.Show;

if nargout == 0
    Show = true;
elseif nargin < 2
    tol = 1e-2;
    % error('chech_massConservation.m: Missing argument of tolerance (pos 2)');
end

[ny, nx, nt] = size(rho);

rho2 = reshape(rho, nx*ny, nt);
sumRho = integralL2(rho2);

negaRho = zeros(size(rho2));
negaRhoInd = rho2 < 0;
negaRho(negaRhoInd) = rho2(negaRhoInd);
sumNegaRho = integralL2(negaRho);

if nargout > 0
    err = norm(cat(1,sumRho-1,sumNegaRho), Inf);
    
    if err > tol
        flag = 0;
    else
        flag = 1;
    end
end

if Show
    disp("Sum of rho: ");
    disp(reshape(sumRho, 1, []));

    disp("Sum of rho(rho<0): ");
    disp(reshape(sumNegaRho, 1, []));
end

end


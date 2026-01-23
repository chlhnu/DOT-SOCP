function var = interpolate(var, model)
%% Interpolates variables (phi and beta) onto a finer grid

% Grid size
nx = model.nx;
nt = model.nt;

% Interpolate
var.phi   = interpolate_phi(var.phi, nx, nt);
var.beta  = interpolate_z(var.beta, nx, nt);

end

%% Interpolate along x/y/t direction
function fR = interpolate_linear(f, d)
    fR = movmean(f, 2, d, 'Endpoints', 'discard');
end

%% Interpolate (t-staggered grid)
function fR = interpolate_tStagger(f)
    % grid size
    [nx, nt] = size(f);
    nxR = 2 * (nx-1) + 1;
    ntR = 2 * nt;
    
    % index
    oddX = 1 : 2 : nxR;
    oddT = 1 : 2 : ntR-1;
    
    evenX = 2 : 2 : nxR-1;
    evenT = 2 : 2 : ntR;

    % interpolate
    fR = zeros(nxR, ntR);
    fR(oddX, oddT) = f;
    fR(oddX, evenT) = f;
    fR(evenX, :) = interpolate_linear(fR(oddX, :), 1);
end

%% Interpolate phi
function phiR = interpolate_phi(phi, nx, nt)
    % grid size
    nxR = 2 * (nx-1) + 1;
    ntR = 2 * (nt-1) + 1;

    % index
    oddX = 1 : 2 : nxR;
    evenX = 2 : 2 : nxR-1;

    oddT = 1 : 2 : ntR;
    evenT = 2 : 2 : ntR-1;
    
    % interpolate
    phiR = zeros(nxR, ntR);
    phiR(oddX , oddT ) = reshape(phi, nx, nt);
    phiR(evenX, oddT ) = interpolate_linear(phiR(oddX, oddT), 1);
    phiR(:    , evenT) = interpolate_linear(phiR(:, oddT), 2);
    
    % reshape
    phiR = reshape(phiR, [], 1);
end

%% Interpolate z
function zR = interpolate_z(z, nx, nt)
    % grid size
    nxR = 2 * (nx-1) + 1;
    ntR = 2 * (nt-1);
    zR = zeros(nxR * ntR, size(z, 2));
    
    for j = 1 : size(z, 2)
        zR(:, j) = reshape( interpolate_tStagger(reshape(z(:, j), nx, nt-1)) , [] , 1);
    end
end


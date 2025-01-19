function var = interpolate(var, model)
%% Interpolates variables (phi and beta) onto a finer grid

% Grid size
nx = model.nx;
ny = model.ny;
nt = model.nt;

% Interpolate
var.phi   = interpolate_phi(var.phi, ny, nx, nt);
var.beta  = interpolate_z(var.beta, ny, nx, nt);

end

%% Interpolate along x/y/t direction
function fR = interpolate_linear(f, d)
    fR = movmean(f, 2, d, 'Endpoints', 'discard');
end

%% Interpolate (t-staggered grid)
function fR = interpolate_tStagger(f)
    % grid size
    [ny, nx, nt] = size(f);
    nyR = 2 * (ny-1) + 1;
    nxR = 2 * (nx-1) + 1;
    ntR = 2 * nt;
    
    % index
    oddY = 1 : 2 : nyR;
    oddX = 1 : 2 : nxR;
    oddT = 1 : 2 : ntR-1;
    
    evenY = 2 : 2 : nyR-1;
    evenX = 2 : 2 : nxR-1;
    evenT = 2 : 2 : ntR;

    % interpolate
    fR = zeros(nyR, nxR, ntR);
    fR(oddY, oddX, oddT) = f;
    fR(oddY, oddX, evenT) = f;
    fR(evenY, oddX, :) = interpolate_linear(fR(oddY, oddX, :), 1);
    fR(:, evenX, :) = interpolate_linear(fR(:, oddX, :), 2);
end

%% Interpolate phi
function phiR = interpolate_phi(phi, ny, nx, nt)
    % grid size
    nyR = 2 * (ny-1) + 1;
    nxR = 2 * (nx-1) + 1;
    ntR = 2 * (nt-1) + 1;

    % index
    oddY = 1 : 2 : nyR;
    evenY = 2 : 2 : nyR-1;

    oddX = 1 : 2 : nxR;
    evenX = 2 : 2 : nxR-1;

    oddT = 1 : 2 : ntR;
    evenT = 2 : 2 : ntR-1;
    
    % interpolate
    phiR = zeros(nyR, nxR, ntR);
    phiR(oddY , oddX , oddT ) = reshape(phi, ny, nx, nt);
    phiR(evenY, oddX , oddT ) = interpolate_linear(phiR(oddY, oddX, oddT), 1);
    phiR(:    , evenX, oddT ) = interpolate_linear(phiR(:   , oddX, oddT), 2);
    phiR(:    , :    , evenT) = interpolate_linear(phiR(:   , :   , oddT), 3);
    
    phiR = reshape(phiR, [], 1);
end

%% Interpolate z
function zR = interpolate_z(z, ny, nx, nt)
    % grid size
    nyR = 2 * (ny-1) + 1;
    nxR = 2 * (nx-1) + 1;
    ntR = 2 * (nt-1);
    zR = zeros(nyR * nxR * ntR, 10);
    
    for j = 1 : 10
        zR(:, j) = reshape( interpolate_tStagger(reshape(z(:, j), ny, nx, nt-1)) , [] , 1);
    end
end


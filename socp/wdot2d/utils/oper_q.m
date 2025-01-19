function diag = oper_q(ny, nx, nt, D, E, weight)
%% Return the elements of the diagonal linear operator I + F^* B^* B F

if ~exist("weight", "var")
    weight = 1;
end

if nargin == 3
    a = repmat(2, [ny, nx, nt-1]);
    b = repmat(2, [ny, nx-1, nt]);
    c = repmat(2, [ny-1, nx, nt]);
    
    b(:, :, [1, end]) = 1;
    c(:, :, [1, end]) = 1;
elseif nargin >= 5
    tmp = (E / D)^2;
    c1  = 2 * tmp;
    c2  = tmp;

    a = repmat(c1, [ny, nx, nt-1]);
    b = repmat(c1, [ny, nx-1, nt]);
    c = repmat(c1, [ny-1, nx, nt]);
    
    b(:, :, [1, end]) = c2;
    c(:, :, [1, end]) = c2;
end

diag = cat(1, a(:), b(:), c(:)) + weight.^2;

end
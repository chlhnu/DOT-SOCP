function diag = oper_q(ny, nx, nt, D, E)
%% Return the elements of diagonal linear operator I + F^* B^* B F
% If the scaling constant (D, E) is provided, 
% return the elements of diagonal linear operator I + (E / D)^2 F^* B^* B F

if nargin == 3 % without scaling
    a = repmat(3, [ny, nx, nt-1]);
    b = repmat(3, [ny, nx-1, nt]);
    c = repmat(3, [ny-1, nx, nt]);
    
    b(:, :, [1, end]) = 2;
    c(:, :, [1, end]) = 2;
elseif nargin == 5 % with scaling
    tmp = (E / D)^2;
    c1  = 1 + 2 * tmp;
    c2  = 1 + tmp;

    a = repmat(c1, [ny, nx, nt-1]);
    b = repmat(c1, [ny, nx-1, nt]);
    c = repmat(c1, [ny-1, nx, nt]);
    
    b(:, :, [1, end]) = c2;
    c(:, :, [1, end]) = c2;
end

diag = cat(1, a(:), b(:), c(:));

end
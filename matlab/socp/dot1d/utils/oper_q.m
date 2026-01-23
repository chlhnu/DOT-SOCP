function diag = oper_q(nx, nt, D, E)
%% The elements of diagonal linear operator I + F^* B^* B F

if nargin == 2 % without scaling
    a = repmat(3, [nx, nt-1]);
    b = repmat(3, [nx-1, nt]);
    b(:, [1, end]) = 2;
elseif nargin == 4 % with scaling
    tmp = (E / D)^2;
    c1  = 1 + 2 * tmp;
    c2  = 1 + tmp;

    a = repmat(c1, [nx, nt-1]);
    b = repmat(c1, [nx-1, nt]);
    b(:, [1, end]) = c2;
else
    error("Input is invalid");
end

diag = cat(1, a(:), b(:));

end
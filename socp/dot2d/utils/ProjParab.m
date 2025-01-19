function q = ProjParab(q)
%% Projection into P

n = size(q, 1);

% helper for cubic equation
extract = @(A)A(:,1);
cubicReal = @(P)real( extract(poly_root(P')') );

% proj
a      = q(:, 1);
b      = q(:, 2:end);
norm_b = vecnorm(b, 2, 2);
lambda = max(cubicReal([ones(n,1), 8 - a, 16 - 8*a, -16*a - 2*norm_b]), 0);
q      = cat(2, a - lambda, b ./ (1 + lambda));

end

function res = normL2(x, h)
%% L2 norm of a column x with stepsize h

res = sqrt(h) * vecnorm(x, 2, 1);

end


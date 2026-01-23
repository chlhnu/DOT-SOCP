function res = FnormL2(x, h)
%% L2 Fnorm of a matrix x with stepsize h

res = sqrt(h) * norm(x, "fro");

end


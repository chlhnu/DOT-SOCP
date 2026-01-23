function res = FnormL2(x, h)
%% L2 Frobenius norm of matrix x with step size h

res = sqrt(h) * norm(x, "fro");

end

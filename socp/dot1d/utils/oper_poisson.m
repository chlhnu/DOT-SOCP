function res = oper_poisson(kernel, rhs)
%% The solution of Poisson equation with Neumann boundary condition

res = reshape( mirt_idctn( mirt_dctn(rhs) ./ kernel ), [], 1);

end

function res = oper_poisson3dim(kernel, rhs)
%% The solution of 3D Poisson equation with Neumann boundary condition

res = reshape( mirt_idctn( mirt_dctn(rhs) ./ kernel ), [], 1);

end

function res = oper_poisson3dim(kernel, rhs)
%% The solution of a 3D Poisson equation with Neumann boundary conditions

res = reshape( mirt_idctn( mirt_dctn(rhs) ./ kernel ), [], 1);

end

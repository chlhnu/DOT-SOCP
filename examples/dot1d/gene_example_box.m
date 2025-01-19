function [rho0, rho1] = gene_example_box(nx)
%% Generate box densities

% Settings
lbox = [0.1, 0.5];
rbox = [0.85, 0.95];

x = reshape(linspace(0,1,nx), nx, 1);
rho0 = double(all(cat(2, x >= lbox(1), x <= lbox(2)), 2));
rho1 = double(all(cat(2, x >= rbox(1), x <= rbox(2)), 2));

end

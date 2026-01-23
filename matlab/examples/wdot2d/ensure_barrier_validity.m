function [rho0, rho1, barrierCenterGrid] = ensure_barrier_validity(rho0, rho1, barrier)
%% Ensure the validity of the barrier in the Weighted DOT problem

[ny, nx] = size(rho0);
[xx, yy] = meshgrid(linspace(0, 1, nx), linspace(0, 1, ny));
barrierCenterGrid = barrier(xx', yy');
filterVal = mean(barrierCenterGrid, 'all');
barrierCenterGrid = (barrierCenterGrid' > filterVal);

rho0(barrierCenterGrid) = 0;
rho1(barrierCenterGrid) = 0;

rho0 = ( (nx * ny / sum(rho0, 'all')) * rho0);
rho1 = ( (nx * ny / sum(rho1, 'all')) * rho1);

end


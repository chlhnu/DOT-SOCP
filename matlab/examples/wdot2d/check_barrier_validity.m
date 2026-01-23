function [barrierCenterGrid, validYes] = check_barrier_validity(rho0, rho1, barrier)
%% Check the validity of the barrier in the Weighted DOT problem

validTol = 1e-4;

[ny, nx] = size(rho0);
[xx, yy] = meshgrid(linspace(0, 1, nx), linspace(0, 1, ny));
barrierCenterGrid = barrier(xx, yy);

if (sum(rho0(barrierCenterGrid)) + sum(rho1(barrierCenterGrid)) > validTol)
    validYes = false;
    error("Invalid (rho0, rho1, barrier)");
else
    validYes = true;
end

end


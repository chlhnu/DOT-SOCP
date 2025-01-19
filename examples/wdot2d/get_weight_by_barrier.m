function weight = get_weight_by_barrier(nx, ny, nt, barrier, barrierWeight)
%% Generate weights based on barrier locations

if ~exist("barrier", "var")
    barrier = @(x,y) false(length(y), length(x));
end

if ~exist("barrerWeight", "var")
    barrierWeight = 1e6;
end

% Grid
hx = 1/(nx-1);
hy = 1/(ny-1);
xStagGrid = linspace(.5*hx, 1-.5*hx, nx-1);
xCentGrid = linspace(0, 1, nx);
yStagGrid = linspace(.5*hy, 1-.5*hy, ny-1);
yCentGrid = linspace(0, 1, ny);

[xx, yy] = meshgrid(xStagGrid, yCentGrid);
mask = (barrier(xx', yy') > 0);
weightX = ones(ny, nx-1);
weightX(mask') = barrierWeight;

[xx, yy] = meshgrid(xCentGrid, yStagGrid);
mask = (barrier(xx', yy') > 0);
weightY = ones(ny-1, nx);
weightY(mask') = barrierWeight;

weightT = ones(ny, nx, nt-1);
weightX = repmat(weightX, [1, 1, nt]);
weightY = repmat(weightY, [1, 1, nt]);
weight = cat(1, weightT(:), weightX(:), weightY(:));

end


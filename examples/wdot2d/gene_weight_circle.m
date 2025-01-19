function weight = gene_weight_circle(nt, nx, ny)
%% Generate weights according to a circular distribution centered at (0.5, 0.5) in a staggered grid

% Circular distribution
[a, b] = deal(.5, .5);
circular = @(xx, yy) sqrt( (xx - a).^2 + (yy - b).^2 );

% Grid
hx = 1/(nx-1);
hy = 1/(ny-1);
xStagGrid = linspace(.5*hx, 1-.5*hx, nx-1);
xCentGrid = linspace(0, 1, nx);
yStagGrid = linspace(.5*hy, 1-.5*hy, ny-1);
yCentGrid = linspace(0, 1, ny);

[xx, yy] = meshgrid(xStagGrid, yCentGrid);
weightX  = circular(xx, yy);
weightX  = weightX * (ny*(nx-1) / sum(weightX, 'all'));

[xx, yy] = meshgrid(xCentGrid, yStagGrid);
weightY  = circular(xx, yy);
weightY  = weightY * (ny*(nx-1) / sum(weightY, 'all'));

weightT = ones(ny, nx, nt-1);
weightX = repmat(weightX, [1, 1, nt]);
weightY = repmat(weightY, [1, 1, nt]);
weight = cat(1, weightT(:), weightX(:), weightY(:));

end

function barrier = gene_barrier_of_example6()
%% Generate labyrinth-shaped obstacles for Example 5.6

barrier = imread('maze.png');
maxVal  = max(barrier, [], 'all');
barrier = (maxVal - barrier) * (255 / maxVal);

[nx, ny] = size(barrier);

[xx, yy] = ndgrid(linspace(0, 1, nx), linspace(0, 1, ny));
barrier = griddedInterpolant(xx, yy, double(barrier), "nearest");

end


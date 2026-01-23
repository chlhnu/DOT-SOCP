function barrier = gene_barrier_of_maze14()
%% Generate a labyrinth-shaped obstacle, which is a labyrinth from [Optimal Transport with Proximal Splitting. SIAM Journal on Imaging Sciences, 2014.]

thisDir = fileparts(mfilename('fullpath'));
barrier = imread(fullfile(thisDir, 'resources', 'maze-14.png'));
maxVal  = max(barrier, [], 'all');
barrier = (maxVal - barrier) * (255 / maxVal);

[nx, ny] = size(barrier);

[xx, yy] = ndgrid(linspace(0, 1, nx), linspace(0, 1, ny));
barrier = griddedInterpolant(xx, yy, double(barrier), "nearest");

end


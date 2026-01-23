function barrier = gene_barrier_of_circle_pillar()
%% Generate a circle-pillar-shaped obstacle

barrier = @(x, y) ( ...
    (x >= 0.2 & x <= 0.25) & ((y >= 0.4 & y <= 1.0)) ... % left pillar
    | (x >= 0.75 & x <= 0.8) & ((y >= 0.0 & y <= 0.6)) ... % right pillar
    | ((x - 0.5).^2 + (y - 0.5).^2 <= 0.15^2) ... % circle
);

end

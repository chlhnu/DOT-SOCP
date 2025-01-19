function [rho0, rho1] = get_example(Problem, nx, ny, lowerBound)
%% Get normalized example

if ~exist('lowerBound', 'var')
    lowerBound = 0;
end

% Generate
if strcmp(Problem, "example1")
    [rho0, rho1] = gene_example1(nx, ny);
elseif strcmp(Problem, "example2")
    [rho0, rho1] = gene_example2(nx, ny);
elseif strcmp(Problem, "example3")
    [rho0, rho1] = gene_example3(nx, ny);
elseif strcmp(Problem, "example4")
    [rho0, rho1] = gene_example4(nx, ny);
elseif strcmp(Problem, "circle")
    [rho0, rho1] = gene_exampleCircle(nx, ny);
elseif strcmp(Problem, "example6")
    [rho0, rho1] = gene_example6(nx, ny);
elseif strcmp(Problem, "maze14")
    [rho0, rho1] = gene_exampleMaze14(nx, ny);
else
    error("Novalid input: 'Problem'");
end

% Add lower bound, normalization
rho0 = ( (nx * ny / sum(rho0, 'all')) * rho0 + lowerBound ) / (1 + lowerBound);
rho1 = ( (nx * ny / sum(rho1, 'all')) * rho1 + lowerBound ) / (1 + lowerBound);

end
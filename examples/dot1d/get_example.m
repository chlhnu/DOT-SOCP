function [rho0, rho1] = get_example(Problem, nx, lowerBound)
%% Get normalized example (1 dimension)

if ~exist('lowerBound', 'var')
    lowerBound = 0;
end

% Generate
if strcmp(Problem, "gaussian")
    [rho0, rho1] = gene_example_gaussian(nx);
elseif strcmp(Problem, "box")
    [rho0, rho1] = gene_example_box(nx);
else
    error("Novalid input: 'Problem'");
end

% Add lower bound, normalization
rho0 = ( (nx / sum(rho0, 'all')) * rho0 + lowerBound ) / (1 + lowerBound);
rho1 = ( (nx / sum(rho1, 'all')) * rho1 + lowerBound ) / (1 + lowerBound);

end
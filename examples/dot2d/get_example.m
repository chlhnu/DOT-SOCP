function [rho0, rho1] = get_example(Problem, nx, ny, lowerBound, varargin)
%% Get normalized example

if ~exist('lowerBound', 'var')
    lowerBound = 0;
end

% If the caller omits lowerBound and starts with name-value pairs, treat the
% first string input as part of varargin instead of lowerBound.
if exist('lowerBound', 'var') && ~isnumeric(lowerBound)
    varargin = [{lowerBound}, varargin];
    lowerBound = 0;
end

parser = inputParser;
addParameter(parser, "DOTmark_type", "ClassicImages");
addParameter(parser, "stitch1_indices", [1, 2, 3, 4]);
addParameter(parser, "stitch2_indices", [5, 6, 7, 8]);
parse(parser, varargin{:});
args = parser.Results;


% Generate
if strcmp(Problem, "example1")
    [rho0, rho1] = gene_example1(nx, ny);
elseif strcmp(Problem, "example2")
    [rho0, rho1] = gene_example2(nx, ny);
elseif strcmp(Problem, "example3")
    [rho0, rho1] = gene_example3(nx, ny);
elseif strcmp(Problem, "example4")
    [rho0, rho1] = gene_example4(nx, ny);
elseif strcmp(Problem, "example5")
    [rho0, rho1] = gene_example5(nx, ny);
elseif strcmp(Problem, "example7")
    [rho0, rho1] = gene_example7(nx, ny);
elseif strcmp(Problem, "circle")
    [rho0, rho1] = gene_exampleCircle(nx, ny);
elseif strcmp(Problem, "DOTmark_4stitch")
    [rho0, rho1] = gene_example_DOTmark_4stitch(nx, ny, args.DOTmark_type, args.stitch1_indices, args.stitch2_indices);
else
    error("Novalid input: 'Problem'");
end

% Add lower bound, normalization
rho0 = ( (nx * ny / sum(rho0, 'all')) * rho0 + lowerBound ) / (1 + lowerBound);
rho1 = ( (nx * ny / sum(rho1, 'all')) * rho1 + lowerBound ) / (1 + lowerBound);

end

function [weight, barrier] = get_weight_from_barrier(Problem, nx, ny, nt)
% Get weight and barrier for predefined WDOT problems
%
%   [weight, barrier] = get_weight_from_barrier(Problem, nx, ny, nt)
%   Returns weight matrix and barrier function handle based on problem type.
%
% INPUTS:
%   Problem - Problem name string:
%             "example8", "example9", "circle-pillar", "maze14"
%   nx, ny  - Spatial grid dimensions
%   nt      - Time grid dimension
%
% OUTPUTS:
%   weight  - Weight matrix (ny, nx, nt)
%   barrier - Barrier function handle

% Generate barrier based on problem type
if strcmp(Problem, "example8")
    barrier = gene_barrier_of_example8();
elseif strcmp(Problem, "example9")
    barrier = gene_barrier_of_example9();
elseif strcmp(Problem, "circle-pillar")
    barrier = gene_barrier_of_circle_pillar();
elseif strcmp(Problem, "maze14")
    barrier = gene_barrier_of_maze14();
else
    error("Invalid problem type for weighted transport. Try: 'example8', 'example9', 'circle-pillar', 'maze14'");
end

% Generate weight from barrier
weight = gene_weight_by_barrier(nx, ny, nt, barrier);

end
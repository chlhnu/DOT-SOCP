function weight = get_weight(Problem, nx, ny, nt)
% Get weight for predefined WDOT problems
%
%   weight = get_weight_from_barrier(Problem, nx, ny, nt)
%
% INPUTS:
%   Problem - Weight name string:
%             "circle", "circleInv"
%   nx, ny  - Spatial grid dimensions
%   nt      - Time grid dimension
%
% OUTPUTS:
%   weight  - Weight matrix (ny, nx, nt)

% Generate barrier based on problem type
if strcmp(Problem, "circle")
    weight = gene_weight_circle(nt, nx, ny);
elseif strcmp(Problem, "circleInv")
    weight = gene_weight_circleInv(nt, nx, ny);
else
    error("Invalid problem type for weighted transport. Try: 'circle', 'circleInv'");
end

end
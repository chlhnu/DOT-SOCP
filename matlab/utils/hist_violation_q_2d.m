function fig = hist_violation_q_2d(q0, bx, by, title, export_yes)
%% Violation of f(q) <= 0

fq = (q0 + .5 * (bx.^2 + by.^2)) ./ (1 + abs(q0) + sqrt(bx.^2 + by.^2));

% Check title
if ~exist('title', 'var') || isempty(title)
    title = "Histogram of f(q) more than 0";
end

% Check export_yes
if ~exist('export_yes', 'var') || isempty(export_yes)
    export_yes = false;
end

fig = hist_positive_value(fq, title, export_yes);

end

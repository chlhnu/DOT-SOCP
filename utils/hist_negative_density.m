function fig = hist_negative_density(rho, title, export_yes)
%% Violation of rho >= 0

% Check title
if ~exist('title', 'var') || isempty(title)
    title = "Histogram of rho less than 0";
end

% Check export_yes
if ~exist('export_yes', 'var') || isempty(export_yes)
    export_yes = false;
end

fig = hist_positive_value(-rho, title, export_yes);

end

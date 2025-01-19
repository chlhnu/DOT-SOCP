function [] = show_residualCurve(resiLists, figName, legendNames, varargin)
%% Plot curve of KKT error

p = inputParser;
addParameter(p, 'xTime', []);
addParameter(p, 'xIteration', []);
addParameter(p, 'lowerBound', NaN);
addParameter(p, 'yLabel', 'Karush-Kuhn-Tucker Errors');
parse(p, varargin{:});
xTime = p.Results.xTime;
xIteration = p.Results.xIteration;
lBound = p.Results.lowerBound;
yLabel = p.Results.yLabel;

num_lgd_items = length(legendNames);
if (num_lgd_items > 0) && (size(resiLists, 2) ~= num_lgd_items)
    error("Error: input data 'resiLists' is not matched with 'legendNames'.");
end

%% Y-axis limit
if isnan(lBound)
    filter = 1e-10;
    lBound = 0.5 * min(resiLists(resiLists > filter), [], 'all');
end

%% Plot
fig = figure('Name', figName);

if ~isempty(xIteration) % x-axis is iteration
    semilogy(xIteration, resiLists, 'LineWidth', 1.5);
    xlabel("Iteration Numbers");
elseif ~isempty(xTime) % x-axis is time
    semilogy(xTime, resiLists, 'LineWidth', 1.5);
    xlabel("Time [seconds]");
else
    semilogy(resiLists, 'LineWidth', 1.5);
end

ylim([lBound Inf]);
ylabel(yLabel);
xlim("tight");

if num_lgd_items >= 4
    lgd = legend(legendNames);
    set(lgd, 'Interpreter', 'latex', 'Location', 'northeastoutside');
    fig_position = [8, 4];
elseif num_lgd_items >= 1
    lgd = legend(legendNames);
    set(lgd, 'Interpreter', 'latex', 'Location', 'best');
    fig_position = [6, 4];
else
    fig_position = [6, 4];
end

set(fig, 'Color', 'w');
set(fig, 'Units', 'inches');
fig.Position(3:4) = fig_position;

ax = gca;
set(ax, 'Fontsize', 12);
set(ax, 'LabelFontSizeMultiplier', 14 / 12);
set(ax, 'TitleFontSizeMultiplier', 16 / 12);
set(ax, 'LineWidth', 1.5);
set(ax, 'FontName', 'Times New Roman');

end

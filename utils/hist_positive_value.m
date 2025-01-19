function [fig, y1_lim, y2_lim] = hist_positive_value(y, title, export_yes, y1_lim, y2_lim)
%% Histogram of y >= 0

y = max(y, 0);

if ~exist('title', 'var') || isempty(title)
    title = "Histogram of value >= 0";
end

if ~exist('export_yes', 'var') || isempty(export_yes)
    export_yes = false;
end

if ~exist('y1_lim', 'var') && ~exist('y2_lim', 'var')
    wo_y_lim = true;
else
    wo_y_lim = false;
end

if (export_yes)
    fig  = figure("Name", title, "Visible", "off");
    path = "results\";
else
    fig = figure("Name", title);
end

%% Histogram
x_left  = -8;
x_right = -1;
levels  = log10(power(10, linspace(x_left, x_right, 50)));
levels_ticks = x_left : 1 : x_right;

log_nega_rho = reshape(log10(y), [], 1);

color_hist = [220,  94,  40] / 255;
color_area = [  0, 114, 189] / 255;
alpha_area = 0.75;

% Y Left
yyaxis left
set(gca, "YColor", color_hist);
histogram(log_nega_rho, levels, 'FaceColor', color_hist, Normalization = "percentage");
ytickformat("percentage");
ylabel("Percentage")

if (~ wo_y_lim)
    ylim([0, y2_lim]);
end

% Y Right
yyaxis right
set(gca, "YColor", color_area);
proportion = histcounts(log_nega_rho, levels, Normalization="percentage");
proportion_cum = cumsum(proportion, 'reverse');
levels = movmean(levels, 2, "Endpoints", "discard");
area(levels, proportion_cum, 'FaceColor', color_area, 'FaceAlpha', alpha_area, 'LineStyle', 'none');
ytickformat("percentage");
ylabel("Cumulative percentage")

if (wo_y_lim)
    y1_lim = proportion_cum(1);
    y2_lim = max(proportion);
else
    ylim([0, y1_lim]);
end

% X
xlabel("Violation");
xticks(levels_ticks);
xticklabels(arrayfun(@(x) sprintf("$$10^{%1d}$$", x), levels_ticks));
set(gca, 'XDir', 'reverse');

% Adjust
fig = adjust_fig(fig);

% Export
if (export_yes)
    exportgraphics(fig, path + title + '.pdf', 'ContentType', 'vector', 'Resolution', 300);
end

end

function fig = adjust_fig(fig)
    axes = findall(fig, 'Type', "Axes");
    set(axes, 'Fontsize', 12);
    set(axes, 'LabelFontSizeMultiplier', 14 / 12);
    set(axes, 'TitleFontSizeMultiplier', 16 / 12);
    set(axes, 'LineWidth', 1.5);
    set(axes, 'FontName', 'Times New Roman');
    arrayfun(@(ax) set(ax.XAxis, 'TickLabelInterpreter', 'latex'), axes);

    set(fig, 'Color', 'w');
    set(fig, 'PaperUnits', 'inches');
    set(fig, 'PaperPosition', [0 0 6 4]);
end


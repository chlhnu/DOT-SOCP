function [] = show_evolution_1d(rho, showFunc, figName)
%% Show evolution of rho
% Supported "showFunc" type:
%   "join" (default), "tile"

if ~exist('showFunc', 'var')
    showFunc = "join";
elseif ~ismember(showFunc, ["tile", "join"])
    error("Invalid input at pos 2");
end

if ~exist('figName', 'var') || isempty(figName)
    figName = "Density evolution";
end

[Nx, Nt] = size(rho);
hx = 1 / Nx;
xx = linspace(hx / 2, 1 - hx / 2, Nx);

y_min = 0;
y_max = 1.1 * max(rho, [], 'all');

%% Display
if strcmp(showFunc, "tile")
    fig = figure("Name", figName, 'Units', 'normalized', 'Position', [0.2, 0.3, 0.6, 0.4]);
    num_tile_horizon = 3;
    num_tile_vertical = 3;
    
    tiledlayout(num_tile_horizon, num_tile_vertical);
    
    for t = round(linspace(1, Nt, num_tile_vertical * num_tile_horizon))
        nexttile;
        plot(xx, rho(:, t), 'LineWidth', 1.5);
        title("$t = " + string(t) + "$", 'interpreter', 'latex');
        ylim([y_min, y_max]);
    end
elseif strcmp(showFunc, "join")
    fig = figure("Name", figName);
    
    num_curves = min(Nt, 10);
    tt = round(linspace(1, Nt, num_curves));

    color_list = custom_turbo(num_curves);
    
    for ind = 1 : num_curves
        plot(xx, rho(:, tt(ind)), 'Color', color_list(ind, :), 'LineWidth', 1.5);
        hold on;
    end
    hold off;
    labels = arrayfun(@(t) sprintf("%.2f", 100 * (t-1) / (Nt-1)), tt);
    lgd = legend(arrayfun(@(label) sprintf("$t=%s(\\%%)$", label), labels, 'UniformOutput', false));
    set(lgd, 'Interpreter', 'latex', 'Location', 'best');
    ylim([y_min, y_max]);
end

adjust_fig(fig);

end

function fig = adjust_fig(fig)
    axes = findall(fig, 'Type', "Axes");
    set(axes, 'Fontsize', 12);
    set(axes, 'LabelFontSizeMultiplier', 14 / 12);
    set(axes, 'TitleFontSizeMultiplier', 16 / 12);
    set(axes, 'LineWidth', 1.5);
    set(axes, 'FontName', 'Times New Roman');
    arrayfun(@(ax) set(ax.XAxis, 'TickLabelInterpreter', 'latex'), axes);
    arrayfun(@(ax) set(ax.YAxis, 'TickLabelInterpreter', 'latex'), axes);

    set(fig, 'Color', 'w');
end

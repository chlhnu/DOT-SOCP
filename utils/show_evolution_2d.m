function [] = show_evolution_2d(rho, showFunc, figName, barrier)
%% Show evolution of rho

if ~exist('showFunc', 'var') || isempty(showFunc)
    showFunc = "imshow";
elseif ~ismember(showFunc, ["imshow", "contourf", "contour", "contour3", "mesh"])
    error("Invalid input at pos 2");
end

if ~exist('figName', 'var') || isempty(figName)
    figName = "Density evolution";
end

if ~exist('barrier', 'var') || isempty(barrier)
    barrier = [];
end

timeDisplay = 3;

% Colormap
colors = turbo;
% cutting_pos = 64;
% colors = colors(cutting_pos + 1 : end, :);
% colors = cat(1, generate_colormap([1, 1, 1], colors(1, :), cutting_pos), colors);

% grid
[ny, nx, nt] = size(rho);
[xx, yy] = meshgrid(linspace(0,1,nx), linspace(0,1,ny));

% barrier
if ~ isempty(barrier)
    if all(size(barrier) == [ny, nx])
        barrier = repmat(barrier, [1, 1, nt]);

        if strcmp(showFunc, "imshow")
            rho(barrier) = Inf;
        elseif strcmp(showFunc, "contour3")
            rho(barrier) = max(rho, [], "all");
        elseif ismember(showFunc, ["contourf", "contour"])
            rho(barrier) = - Inf;
        elseif strcmp(showFunc, "mesh")
            error("Type of mesh is invalid if there is barrier");
        end

        colors = turbo;
    else
        error("Argument at posistion 5 is invalid");
    end
end

%% Display type
maxVal     = max(rho, [], "all");

if strcmp(showFunc, "imshow")
    fixAxis    = @() 1;
    disPlayRho = @(t) imshow(rho(:, :, t), [0, maxVal], 'InitialMagnification', 512);
elseif strcmp(showFunc, "contourf")
    rho = rho * (255 / maxVal);
    fixAxis = @() colormap(colors);
    if ~ isempty(barrier)
        levels = [-10, exp(linspace(0, log(255), 128))];
    else
        levels = exp(linspace(0, log(255), 128));
    end
    disPlayRho = @(t) contourf(xx, yy, rho(:, :, t), levels, 'LineColor', 'none');
elseif strcmp(showFunc, "contour")
    rho = rho * (255 / maxVal);
    fixAxis    = @() colormap(colors);
    if ~ isempty(barrier)
        levels  = [-10, exp(linspace(0, log(255), 30))];
    else
        levels  = exp(linspace(0, log(255), 30));
    end
    disPlayRho = @(t) contour(xx, yy, rho(:, :, t), levels);
elseif strcmp(showFunc, "contour3")
    if ~ isempty(barrier)
        fixAxis    = @() view(gca, [-35.1 86.0544759469207]); axis([0, 1, 0, 1, 0, maxVal]);
    else
        fixAxis    = @() axis([0, 1, 0, 1, 0, maxVal]);
    end
    disPlayRho = @(t) contour3(xx, yy, rho(:, :, t), 30);
elseif strcmp(showFunc, "mesh")
    fixAxis    = @() axis([0, 1, 0, 1, 0, maxVal]);
    disPlayRho = @(t) mesh(xx, yy, rho(:, :, t), 'FaceColor', 'flat');
end

%% Display
fig = figure("Name", figName);

for t = 1 : nt
    disPlayRho(t);
    fixAxis();
    adjust_fig(fig);
    pause(timeDisplay / nt);
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
    arrayfun(@(ax) set(ax.YAxis, 'TickLabelInterpreter', 'latex'), axes);

    set(fig, 'Color', 'w');
    set(fig, 'PaperUnits', 'inches');
    set(fig, 'PaperPosition', [0 0 6 4]);
end


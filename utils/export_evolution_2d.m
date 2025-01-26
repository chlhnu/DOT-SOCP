function [ ] = export_evolution_2d(rho, path_to_exported_file, num_frame, showFunc, barrier)
%% Export evolution of rho

if ~exist("num_frame", "var") || isempty(num_frame)
    num_frame = 6;
end

if ~exist('showFunc', 'var') || isempty(showFunc)
    showFunc = "imshow";
elseif ~ismember(showFunc, ["imshow", "contourf", "contour", "contour3", "mesh"])
    error("Invalid input at pos 4");
end

if ~exist('barrier', 'var') || isempty(barrier)
    barrier = [];
end

%% Exported file

default_exported_path = "results/";

% Check <path_to_exported_file>
if ~exist("path_to_exported_file", "var")
    timeStr = string(datetime('now','TimeZone','local','Format','yyyy-MM-dd-HH-mm-ss'));
    path_to_exported_file = fullfile(default_exported_path, timeStr + ".pdf");
end

[exported_path, exported_filename, ext] = fileparts(path_to_exported_file);

% Creat directory
[status, ~] = mkdir(exported_path);
if status == 0
    error("Failed to create the directory for the exported file.")
end

valid_image_ext = [".pdf", ".png", ".jpg"];
valid_video_ext = [".mp4", ".avi"];

if ismember(ext, valid_image_ext)
    to_export_img = true;
    to_export_vdo = false;
elseif ismember(ext, valid_video_ext)
    [~, ~, nt] = size(rho);
    num_frame = nt;

    to_export_img = false;
    to_export_vdo = true;
else
    error("The filename extension of input at position 2 is invalid");
end


%% Plot
% Colormap
colors = custom_turbo();
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

        colors = custom_turbo();
    else
        error("Argument at posistion 5 is invalid");
    end
end

%% Display type
maxVal = max(rho, [], "all");

if strcmp(showFunc, "imshow")
    fixAxis    = @() 1;
    setGraphic = @() setGraphic_WOAxis();
    rho2       = maxVal - rho;
    disPlayRho = @(t) imshow(rho2(:, :, t), [0, maxVal]);
elseif strcmp(showFunc, "contourf")
    rho = rho * (255 / maxVal);
    fixAxis    = @() colormap(colors);
    setGraphic = @() setGraphic_WO2dimAxis();
    if ~ isempty(barrier)
        levels = [-10, exp(linspace(0, log(255), 128))];
    else
        levels = exp(linspace(0, log(255), 128));
    end
    disPlayRho = @(t) contourf(xx, yy, rho(:, :, t), levels, 'LineColor', 'none');
elseif strcmp(showFunc, "contour")
    rho = rho * (255 / maxVal);
    fixAxis    = @() colormap(colors);
    setGraphic = @() setGraphic_W2dimAxis();
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
    setGraphic = @() setGraphic_W3dimAxis();
    disPlayRho = @(t) contour3(xx, yy, rho(:, :, t), 30);
elseif strcmp(showFunc, "mesh")
    fixAxis    = @() axis([0, 1, 0, 1, 0, maxVal]);
    setGraphic = @() setGraphic_W3dimAxis();
    disPlayRho = @(t) mesh(xx, yy, rho(:, :, t), 'FaceColor', 'flat');
end

%% Export graphics to path

ind_frame = round(linspace(1, nt, num_frame));

height = 800; % unit: pixels
width  = 800; % unit: pixels
vdo_width = 1200; % unit: pixels

fig = figure('Color', 'w', 'unit', 'pixels', 'position', [100, 100, width, height]);

if to_export_vdo % Export 4:3 video
    vdo_container = zeros(height, vdo_width, 3, 'uint8');
    
    left_margin = floor((vdo_width - width) / 2);
    vdo_content_idx = zeros(size(vdo_container));
    vdo_content_idx(1 : height, left_margin+1 : left_margin+width, :) = 1;
    vdo_content_idx = logical(vdo_content_idx);
end

if to_export_img
    if num_frame <= 20
        flops = 2;
    elseif num_frame <= 200
        flops = 3;
    else
        flops = 4;
    end
    
    t_frame = round((ind_frame - 1) ./ (nt - 1), flops);
    
    % Ignore the frame of "t=0"
    % num_frame = num_frame - 1;
    % ind_frame = ind_frame(2:end);
    % t_frame   = t_frame(2:end);
elseif to_export_vdo
    if strcmp(ext, ".mp4")
        vdo = VideoWriter(fullfile(exported_path, exported_filename + ext), "MPEG-4");
        vdo.Quality = 100;
    elseif strcmp(ext, ".avi")
        vdo = VideoWriter(fullfile(exported_path, exported_filename + ext), "Motion JPEG AVI");
        vdo.Quality = 80;
    else
        error("An exception occurs");
    end
    total_time = 5; % unit: secondes
    vdo.FrameRate = num_frame / total_time;
    open(vdo);
else
    error("An exception occurs");
end

% Plot
for k = 1 : num_frame
    clf(fig);
    
    % Plot one frame
    ind = ind_frame(k);
    disPlayRho(ind);
    
    % if strcmp(showFunc, "contourf") && ind == nt
    %     contourf(rho(:, :, ind), exp(linspace(0, log(max(rho, [], "all")), 128)), 'LineColor', colors(end, :), 'LineWidth', .1);
    % else
    %     disPlayRho(ind);
    % end

    fixAxis();
    setGraphic();

    % Export
    if to_export_img
        name2 = num2str(t_frame(k), sprintf("%%.%df", flops));
        exportgraphics(fig, fullfile(exported_path, exported_filename + "-t=" + name2 + ext), "Resolution", 600);
    elseif to_export_vdo
        setGraphic_video();
        frame = getframe(fig);
        vdo_container(vdo_content_idx) = imresize(frame.cdata, [width, height]);
        writeVideo(vdo, vdo_container);
    else
        error("An exception occurs");
    end
end

% Delete
delete(fig);
if to_export_vdo
    close(vdo);
end

end

%% Set graphic

function [] = setGraphic_video()
    % set(gcf, 'unit', 'pixels', 'position', [100, 100, 1200, 800]);
    set(gcf, 'Renderer', 'opengl');
    % axis tight;
    % axis padded;
end

function [] = setGraphic_WOAxis()
    leftMargin = 0.05; % Unit: percent
    bottomMargin = 0;  % Unit: percent
    topMargin    = 0; % Unit: percent

    set(gca, 'Position', [leftMargin, bottomMargin, 1-2*leftMargin, 1-bottomMargin-topMargin]);
end

function [] = setGraphic_W2dimAxis()
    leftMargin = 0.15; % Unit: percent
    rightMargin = 0.05; % Unit: percent
    bottomMargin = 0.1;  % Unit: percent
    topMargin    = 0.05; % Unit: percent

    ax = gca;
    set(ax, 'Position', [leftMargin, bottomMargin, 1-leftMargin-rightMargin, 1-bottomMargin-topMargin]);
    ax.XGrid = false;
    ax.YGrid = false;
    ax.ZGrid = false;
    ax.LineWidth = 1;
    ax.FontName = 'Times New Roman';
    ax.FontSize = 12;
    ax.TickLabelInterpreter = 'latex';
end

function [] = setGraphic_WO2dimAxis()
    leftMargin = 0.00; % Unit: percent
    rightMargin = 0.00; % Unit: percent
    bottomMargin = 0.00;  % Unit: percent
    topMargin    = 0.00; % Unit: percent

    box on;
    ax = gca;
    set(ax, 'Position', [leftMargin, bottomMargin, 1-leftMargin-rightMargin, 1-bottomMargin-topMargin]);
    set(ax, 'xtick', [], 'ytick', []);
    ax.XGrid = false;
    ax.YGrid = false;
    ax.LineWidth = 2;
    ax.FontName = 'Times New Roman';
    ax.FontSize = 12;
    ax.TickLabelInterpreter = 'latex';
end

function [] = setGraphic_W3dimAxis()
    leftMargin = 0.00; % Unit: percent
    rightMargin = 0.00; % Unit: percent
    bottomMargin = 0.00;  % Unit: percent
    topMargin    = 0.00; % Unit: percent

    axis off;
    ax = gca;
    set(ax, 'Position', [leftMargin, bottomMargin, 1-leftMargin-rightMargin, 1-bottomMargin-topMargin]);
    set(ax, 'xtick', [], 'ytick', [], 'ztick', []);
    ax.XGrid = false;
    ax.YGrid = false;
    ax.ZGrid = false;
    ax.LineWidth = 1;
    ax.FontName = 'Times New Roman';
    ax.FontSize = 12;
    ax.TickLabelInterpreter = 'latex';
end

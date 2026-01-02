function [rho0, rho1] = gene_example_DOTmark_4stitch(nx, ny, type, stitch1_indices, stitch2_indices)
%% Example 5.6: stitch of images from DOTmark

% Check `type` input
if ~exist("type", "var")
    type = lower("ClassicImages");
else
    type = lower(type);
    if ~ismember(type, ["classicimages", "shapes"])
        error("type (pos 3) must be 'ClassicImages' or 'Shapes'");
    end
end

% Check `stitch1_indices` input
if ~exist("stitch1_indices", "var")
    stitch1_indices = [1, 2, 3, 4];
else
    if length(stitch1_indices) ~= 4
        error("stitch1_indices (pos 4) must be a vector of length 4");
    end
end

% Check `stitch2_indices` input
if ~exist("stitch2_indices", "var")
    stitch2_indices = [5, 6, 7, 8];
else
    if length(stitch2_indices) ~= 4
        error("stitch2_indices (pos 5) must be a vector of length 4");
    end
end

this_dir = fileparts(mfilename("fullpath"));
if type == "classicimages"
    imgs_dir = fullfile(this_dir, "resources", "DOTmark", "ClassicImages");
elseif type == "shapes"
    imgs_dir = fullfile(this_dir, "resources", "DOTmark", "Shapes");
else
    error("Unexpected type");
end

stitch1_paths = {...
    fullfile(imgs_dir, sprintf("%d.png", stitch1_indices(1))), ...
    fullfile(imgs_dir, sprintf("%d.png", stitch1_indices(2))), ...
    fullfile(imgs_dir, sprintf("%d.png", stitch1_indices(3))), ...
    fullfile(imgs_dir, sprintf("%d.png", stitch1_indices(4))) ...
};

stitch2_paths = {...
    fullfile(imgs_dir, sprintf("%d.png", stitch2_indices(1))), ...
    fullfile(imgs_dir, sprintf("%d.png", stitch2_indices(2))), ...
    fullfile(imgs_dir, sprintf("%d.png", stitch2_indices(3))), ...
    fullfile(imgs_dir, sprintf("%d.png", stitch2_indices(4))) ...
};

rho0 = gene_big_image(stitch1_paths, nx, ny);
rho1 = gene_big_image(stitch2_paths, nx, ny);

if size(rho0, 1) ~= ny || size(rho0, 2) ~= nx
    rho0 = imresize(rho0, [ny, nx]);
end
if size(rho1, 1) ~= ny || size(rho1, 2) ~= nx
    rho1 = imresize(rho1, [ny, nx]);
end

rho0 = double(rho0);
rho1 = double(rho1);

end

function rho = gene_big_image(paths, total_h, total_w)
    if length(paths) ~= 4
     error('Please imput 4 images');
    end
    
    sub_h = floor(total_h / 2);
    sub_w = floor(total_w / 2);
    
    imgs = cell(1,4);
    for i = 1:4
        img = imread(paths{i});
        if size(img, 3) == 3
            img = rgb2gray(img);
        end
        img = imresize(img, [sub_h, sub_w]);
        imgs{i} = double(img);
    end

    rho = [imgs{1}, imgs{2}; imgs{3}, imgs{4}];
end

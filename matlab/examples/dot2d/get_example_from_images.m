function [rho0, rho1] = get_example_from_images(path_to_img1, path_to_img2, nx, ny, varargin)
%% Get normalized example based on two images

if any(strcmpi(varargin, 'ReverseColor'))
    is_reverse_color = true;
else
    is_reverse_color = false;
end

% Read images
rho0 = imresize(imread(path_to_img1), [ny-1, nx-1]);
rho1 = imresize(imread(path_to_img2), [ny-1, nx-1]);

if size(rho0, 3) == 3
    rho0 = rgb2gray(rho0);
end

if size(rho1, 3) == 3
    rho1 = rgb2gray(rho1);
end

if is_reverse_color
    rho0 = double(255 - rho0) / 255;
    rho1 = double(255 - rho1) / 255;
end

rho0 = custom_padarray(rho0, [1, 1], 0, 'post');
rho1 = custom_padarray(rho1, [1, 1], 0, 'post');

% Normalization
rho0 = (nx * ny / sum(rho0, 'all')) * double(rho0);
rho1 = (nx * ny / sum(rho1, 'all')) * double(rho1);

end
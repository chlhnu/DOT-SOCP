function [rho0, rho1] = gene_example5(nx, ny)
%% Example 5.5

% read img
rho0 = imresize(imread('centaur.bmp'), [ny-1, nx-1]);
rho1 = imresize(imread('man.bmp'), [ny-1, nx-1]);

% color inversion
rho0 = double(255 - rho0) / 255;
rho1 = double(255 - rho1) / 255;

% padding zero
rho0 = custom_padarray(rho0, [1, 1], 0, 'post');
rho1 = custom_padarray(rho1, [1, 1], 0, 'post');

end


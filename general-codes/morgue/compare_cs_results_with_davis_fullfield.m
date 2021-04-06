clear
close all
clc

restoredefaultpath;
addpath(genpath('/scratch/shannon/c/aether/Projects/BOS/general-codes/matlab-codes/'));
addpath ../prana-uncertainty-average-dc-new-im-2-save-deform-cs/
% addpath ../prana/
%% image settings

% window resolution
window_resolution = 32;

% window size
window_size = 2 * window_resolution;
dataset_name = 'PivChal03B';

%% plot settings

% range of displacements to be displayed in the contour plots
displacement_color_min = 0;
displacement_color_max = 15;
displacement_contour_levels = linspace(displacement_color_min, displacement_color_max, 100);

% range of uncertainties to be displayed (pix.)
uncertainty_color_min = 0;
uncertainty_color_max = 0.1;
uncertainty_contour_levels = linspace(uncertainty_color_min, uncertainty_color_max, 100);
max_bin = 0.08;

%% load davis results

% load a displacement field from davis
davis_results_filename = fullfile('/scratch/shannon/a/aether/Projects/PIV_PTV_Uncertainty_Quantification/Planar_uncertainty_work/Results/07_30_2018/', ['WS' num2str(window_resolution)], '/SOC_gradient/matfiles/', [dataset_name '.mat']);
% results_filename = fullfile('../results/prana-test/PIV_pass4_001.mat');
results = load(davis_results_filename);

snapshot_index = 1;

% extract davis results
Xd = results.Xd;
Yd = results.Yd;

sigma_cs_x_d = results.UCSx(:, :, snapshot_index);
sigma_cs_y_d = results.UCSy(:, :, snapshot_index);


%% load prana results

% results directory for prana
% prana_results_directory = '/scratch/shannon/a/aether/Projects/PIV_PTV_Uncertainty_Quantification/uncertainty-combination/analysis/results/prana-test-cs-buffer=12/';
prana_results_directory = '/scratch/shannon/a/aether/Projects/PIV_PTV_Uncertainty_Quantification/uncertainty-combination/analysis/results/experiment-new/PivChal03B/WS1/vectors/';

% load vector field
% prana_results = load(fullfile(prana_results_directory, 'PIV_pass4_001.mat'));
prana_results = load(fullfile(prana_results_directory, 'PIV_pass4_0001.mat'));

% extract displacement field
Xp = prana_results.X;
Yp = prana_results.Y;

%% display results as contourf plots

figure
subplot(2, 1, 1)
contourf(Xp, Yp, sqrt(prana_results.U.^2 + prana_results.V.^2), displacement_contour_levels, 'linestyle', 'none'), colormap(flipud(gray))
colorbar;
annotate_image(gcf, gca);
title('Prana')

subplot(2, 1, 2)
contourf(Xd, Yd, sqrt(results.Ud(:, :, 1).^2 + results.Vd(:, :, 1).^2), displacement_contour_levels, 'linestyle', 'none'), colormap(flipud(gray))
colorbar;
annotate_image(gcf, gca);
title('Davis')


%% calculate uncertainty

% load deformed images
% im1 = imread(fullfile(prana_results_directory, 'imDeform', 'PIV_pass4_im1d_001.tif'));
% im2 = imread(fullfile(prana_results_directory, 'imDeform', 'PIV_pass4_im2d_001.tif'));
im1 = imread(fullfile(prana_results_directory, 'imDeform', 'PIV_pass4_im1d_0001.tif'));
im2 = imread(fullfile(prana_results_directory, 'imDeform', 'PIV_pass4_im2d_0001.tif'));

% extract image size
[image_height, image_width] = size(im1);

% covert to double
im1 = double(im1);
im2 = double(im2);

% flip images
im1 = flipud(im1);
im2 = flipud(im2);

% deform images based on velocity field
im1def = im1;
im2def = im2;

% calculate uncertainty using cs
% [sigma_cs_x_temp, sigma_cs_y_temp] = correlation_statistics(flipud(im1def), flipud(im2def), [window_size, window_size], [window_resolution, window_resolution; window_resolution, window_resolution], Xd(:), Yd(:), 2, 2); 
[sigma_cs_x_temp, sigma_cs_y_temp] = correlation_statistics(im1def, im2def, [window_size, window_size], [window_resolution, window_resolution; window_resolution, window_resolution], Xp(prana_results.Eval >= 0), Yp(prana_results.Eval >= 0), 2, 2); 

% reshape array to matrix form
sigma_cs_x_p = abs(reshape(sigma_cs_x_temp, size(Xp, 1), size(Xp, 2)));
sigma_cs_y_p = abs(reshape(sigma_cs_y_temp, size(Xp, 1), size(Xp, 2)));

%% display results as contourf plots

figure
subplot(2, 1, 1)
contourf(Xp, Yp, sigma_cs_x_p, uncertainty_contour_levels, 'linestyle', 'none'), colormap(flipud(gray))
colorbar;
annotate_image(gcf, gca);
title('Prana, X')

subplot(2, 1, 2)
contourf(Xd, Yd, sigma_cs_x_d, uncertainty_contour_levels, 'linestyle', 'none'), colormap(flipud(gray))
colorbar;
annotate_image(gcf, gca);
title('Davis, X')

figure
subplot(2, 1, 1)
contourf(Xp, Yp, sigma_cs_y_p, uncertainty_contour_levels, 'linestyle', 'none'), colormap(flipud(gray))
colorbar;
caxis([0 0.1])
annotate_image(gcf, gca);
title('Prana, Y')

subplot(2, 1, 2)
contourf(Xd, Yd, sigma_cs_y_d, uncertainty_contour_levels, 'linestyle', 'none'), colormap(flipud(gray))
colorbar;
annotate_image(gcf, gca);
title('Davis, Y')

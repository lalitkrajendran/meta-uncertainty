clear
close all
clc

addpath(genpath('/scratch/shannon/c/aether/Projects/BOS/general-codes/matlab-codes/'));
addpath('prana-uncertainty-average-dc-new-im-2/');

%% image settings

% window resolution
window_resolution = 64;

% window size
window_size = 2 * window_resolution;
dataset_name = 'PivChal03B';

% load a displacement field from davis
results_filename = fullfile('/scratch/shannon/a/aether/Projects/PIV_PTV_Uncertainty_Quantification/Planar_uncertainty_work/Results/07_30_2018/WS64/SOC_gradient/matfiles/', [dataset_name '.mat']);
results = load(results_filename);

snapshot_index = 1;

%% plot settings
% range of displacements to be displayed in the contour plots
displacement_color_min = 0;
displacement_color_max = 15;
displacement_contour_levels = linspace(displacement_color_min, displacement_color_max, 100);

% maximum error for a measurement to be considered valid
max_error_threshold = 0.1;

% range of uncertainties to be displayed (pix.)
uncertainty_color_min = 0;
uncertainty_color_max = 0.1;
uncertainty_contour_levels = linspace(uncertainty_color_min, uncertainty_color_max, 100);
max_bin = 0.08;

%% uncertainty settings

% use gaussian filter for CS? (True/false)
gaussian_filtering = true;

% size of the Gaussian filter (pix.)
gaussian_filter_width_array = [0, 2, 4, 8, 16, 32, 64];
num_filter_width = length(gaussian_filter_width_array);

% uncertainty flags for the various methods (0 = False, 1 = True)
uncertainty_flags.ppruncertainty = 0;
uncertainty_flags.miuncertainty = 0;
uncertainty_flags.mcuncertainty = 1;

%% load image pair
im1 = imread('/scratch/shannon/a/aether/Projects/PIV_PTV_Uncertainty_Quantification/Planar_uncertainty_work/Images/PivChal03B/B001.tif');
im2 = imread('/scratch/shannon/a/aether/Projects/PIV_PTV_Uncertainty_Quantification/Planar_uncertainty_work/Images/PivChal03B/B002.tif');

% covert to double
im1 = double(im1);
im2 = double(im2);

% flip images
im1 = flipud(im1);
im2 = flipud(im2);

[image_height, image_width] = size(im1);

%% extract displacement field for this snapshot
X = results.Xd;
Y = results.Yd;
U = results.Ud(:, :, snapshot_index);
V = results.Vd(:, :, snapshot_index);                

%% deform images based on velocity field

% create pixel grid
[XI,YI] = IMgrid([image_height, image_width], [0 0]);
XI = XI - 0.5;
YI = YI - 0.5;

% interpolate displacements onto pixel grid
UI = VFinterp(X, Y, U, XI, YI, 3);
VI = VFinterp(X, Y, V, XI, YI, 3);

% translate pixel locations, but
% since sincBlackmanInterp2 assumes
% coordinate system is pixel-centered, we need
% to convert back to index-coordinates for the
% deform.
%HOWEVER: we need to use these shifted
%coordinates later in vector coordinates, so
%move the -0.5 pixel correction to the call to
%the sincBlackmanInterp2 function
XD1 = XI - UI/2;
YD1 = YI - VI/2;
XD2 = XI + UI/2;
YD2 = YI + VI/2;

% im1d = interp2(im1, XD1+0.5, YD1+0.5, 'cubic',0);
% im2d = interp2(im2, XD2+0.5, YD2+0.5, 'cubic',0);

% im1d = whittaker_blackman(im1, XD1 + 0.5, YD1 + 0.5, 6, 1);
% im2d = whittaker_blackman(im2, XD1 + 0.5, YD1 + 0.5, 6, 1);
im1d = whittaker_blackman(im1, XD1, YD1, 6, 1);
im2d = whittaker_blackman(im2, XD1, YD1, 6, 1);

%% select grid points for analysis
r = 20; %[20, 60];
c = 51; %[51, 51];

%% identify window extents
xmin = X(r, c) - window_size/2 + 1;
xmax = X(r, c) + window_size/2;
ymin = Y(r, c) - window_size/2 + 1;
ymax = Y(r, c) + window_size/2;

%% extract image windows

% find the image windows 
zone1 = im1d( max([1 ymin]):min([image_height ymax]),max([1 xmin]):min([image_width xmax]));
zone2 = im2d( max([1 ymin]):min([image_height ymax]),max([1 xmin]):min([image_width xmax]));

if size(zone1,1)~=window_size || size(zone1,2)~=window_size
    w1 = zeros(window_size,window_size);
    w1( 1+max([0 1-ymin]):window_size-max([0 ymax-image_height]), 1+max([0 1-xmin]):window_size-max([0 xmax-image_width]) ) = zone1;
    zone1 = w1;
end
if size(zone2,1)~=window_size || size(zone2,2)~=window_size
    w2 = zeros(window_size,window_size);
    w2( 1+max([0 1-ymin]):window_size-max([0 ymax-image_height]),1+max([0 1-xmin]):window_size-max([0 xmax-image_width]) ) = zone2;
    zone2 = w2;
end

%% apply windowing function to interrogation window

% create window masking filter
sfilt1 = windowmask([window_size window_size], [window_resolution window_resolution]);
sfilt2 = windowmask([window_size window_size], [window_resolution window_resolution]);

% apply the image spatial filter
im1_sub = zone1 .* sfilt1;
im2_sub = zone2 .* sfilt2;

%% calculate uncertainty using IM
[sigma_im_x, sigma_im_y, ~, ~] = original_particle_disparity_no_deform(im1_sub, im2_sub, X(r, c), Y(r, c), U(r, c), V(r, c), image_height);

%% calculate uncertainty using MC
[~,~,u_sub,v_sub,~,~,~,uncertainty2D,~] = PIVwindowed(im1_sub,im2_sub,'SCC',[window_size, window_size], [window_resolution, window_resolution; window_resolution, window_resolution], 0, [2.8, 2.8], 1, 3, 0, 0, 0, X(r, c), Y(r, c), uncertainty_flags, 1, 0, 0);

sigma_mc_x = sqrt(uncertainty2D.biasx.^2 + (uncertainty2D.Ixx.^2)./uncertainty2D.Neff);
sigma_mc_y = sqrt(uncertainty2D.biasy.^2 + (uncertainty2D.Iyy.^2)./uncertainty2D.Neff);

U_new = u_sub;
V_new = v_sub;
%% calculate cs uncertainty
sigma_cs_x = nans(1, num_filter_width);
sigma_cs_y = nans(1, num_filter_width);

for filter_width_index = 1:num_filter_width    
    [sigma_cs_x(filter_width_index), sigma_cs_y(filter_width_index)] = correlation_statistics(im1_sub, im2_sub, [128, 128], [128, 128], 64, 64, true, gaussian_filter_width_array(filter_width_index));
end

% [x, y] = correlation_statistics_original(im1_sub, im2_sub, [128, 128], [128, 128], 64, 64)
%% display results
fprintf('Davis: %.3f\n', results.UCSx(r, c, snapshot_index));
fprintf('Prana: %.3f\n', sigma_cs_x);

%%
colors = lines(2);
figure
plot(gaussian_filter_width_array, sigma_cs_x, '*', 'color', colors(1,:))
hold on
plot([0, 64], results.UCSx(r, c, snapshot_index) * [1, 1], 'color', colors(1,:))
plot(gaussian_filter_width_array, sigma_cs_y, '*', 'color', colors(2,:))
plot([0, 64], results.UCSy(r, c, snapshot_index) * [1, 1], 'color', colors(2,:))
plot([0, 64], results.UCSx(r, c, snapshot_index) * [1, 1], 'color', colors(1,:))

xlabel('Filter Width')
ylabel('\sigma_{CS} (pix.)')
legend('x - Prana', 'x - Davis', 'y - Prana', 'y- Davis', 'location', 'eastoutside')
% This script calculates the pdf of uncertainties estimated by three
% different algorithms for experimental images by jack-knife resampling.

%%
clear
close all
clc

restoredefaultpath;
addpath(genpath('/scratch/shannon/c/aether/Projects/BOS/general-codes/matlab-codes/'));
addpath('../prana-uncertainty-average-dc-new-im-2-save-deform-cs/');
addpath('../general-codes/')
setup_default_settings;
% dbstop if error

% ====================================
%% read/write settings
% ====================================
% top level read directory
top_read_directory = '/scratch/shannon/a/aether/Projects/PIV_PTV_Uncertainty_Quantification/Stereo_uncertainty_work/';
% top level write directory
top_write_directory = '/scratch/shannon/a/aether/Projects/PIV_PTV_Uncertainty_Quantification/uncertainty-combination/analysis/results/stereo-dataset/cam13/With_selfcal/';
mkdir_c(top_write_directory);
% Load 2d job file 
job_settings = load(fullfile(top_write_directory, 'cam13_VR_prana_fulljob_withselfcal.mat'));
% extract calibration job file
caljobfile = job_settings.caljobfile;
% extract 2d job file
planarjob = job_settings.planarjob;

% ====================================
%% processing settings
% ====================================
% camera numbers
camera_numbers = [1, 3];
% stereo reconstruction type
rectype = 'Willert';
% pass number
pass_number = 4;
% array of uncertainty methods
uncertainty_method_array = {'IM'; 'MC'; 'CS'};
num_uncertainty_methods = numel(uncertainty_method_array);
% start frame
start_frame = 25;

% ====================================
%% resampling settings
% ====================================
% number of particles to remove
% num_particles_remove_array = [5, 10, 15];
percentage_particles_remove = 0.25;
% number of resampling trials
num_resampling_trials = 1e2;
% uncertainty flags for the various methods (0 = False, 1 = True)
uncertainty_flags.ppruncertainty = 0;
uncertainty_flags.miuncertainty = 0;
uncertainty_flags.mcuncertainty = 1;

% ====================================
%% particle identification settings
% ====================================
% intensity threshold for particle identification
intensity_threshold = 10;
% particle diameter (pix.)
d_p = 3;
% particle sizing settings
sizeprops = struct;
sizeprops.method = 'IWC';
sizeprops.p_area = 2; % this will be modified for each image
sizeprops.sigma = 4;
sizeprops.errors = 0;
sizeprops.save_dir = [];
% display identification results? 
display_id_results = 0;

% ====================================
%% load 2d, 2c vector fields
% ====================================
fprintf('Loading all datasets into memory\n');

results_2c_all = cell(1, 2);
jobfile_all = cell(1, 2);
files_im1 = cell(1, 2);
files_im2 = cell(1, 2);

% loop through cameras
for camera_index = 1:2
    fprintf('Camera: %d\n', camera_index);
    % ====================================
    %% Load data
    % ====================================
    % results directory for the current camera
    current_results_directory = fullfile(top_write_directory, rectype, ['Camera', num2str(camera_index), filesep], 'vectors');
    % load results for vectors, errors and uncertainties
    [results_2c_all{camera_index}, num_snapshots] = load_directory_data(current_results_directory, ['VR*pass' num2str(pass_number, '%d') '*.mat']);
    
    % ====================================
    %% load listing of deformed images
    % ====================================
    
    % directory containing deformed images for current data set
    deformed_images_directory = fullfile(current_results_directory, 'imDeform');
    
    % get list of im1 files in the directory
    files_im1{camera_index} = get_directory_listing(deformed_images_directory, 'VR*im1d*.tif');
    % get list of im2 files in the directory
    files_im2{camera_index} = get_directory_listing(deformed_images_directory, 'VR*im2d*.tif');

end

% ====================================
%% load errors and uncertainties
% ====================================
% directory containing stereo results
stereo_results_directory = fullfile(top_write_directory, rectype, 'Camera1Camera2_3Cfields',filesep);

% load errors
fprintf('Loading all errors into memory\n');
filename = fullfile(stereo_results_directory, 'errors.mat');
errors_all = load(filename);

% load uncertainties
fprintf('Loading all uncertainties into memory\n');
filename = fullfile(stereo_results_directory, 'uncertainties.mat');
uncertainties_all = load(filename);

% ====================================
%% run resampling
% ====================================

for snapshot_index = 1 %:num_snapshots
    % loop through cameras
    for camera_index = 1:2 %1:2
        % load camera 1 result
        current_result = results_2c_all{camera_index}(snapshot_index);        
        % extract number of grid point
        [num_rows, num_cols] = size(current_result.X);
        num_grid_points = num_rows * num_cols;

        % extract job file
        jobfile = job_settings.(['job' num2str(camera_index)]);

        % ================================================
        % load deformed images
        % ================================================
        im1 = imread(fullfile(files_im1{camera_index}(snapshot_index).folder, files_im1{camera_index}(snapshot_index).name));
        im1d = double(flipud(im1));
        
        im2 = imread(fullfile(files_im2{camera_index}(snapshot_index).folder, files_im2{camera_index}(snapshot_index).name));
        im2d = double(flipud(im2));
        
        % extract image height and width
        [image_height, image_width] = size(im1d);
        
        % initialize variables to store resampling results
        unc_resampling = cell(1, num_grid_points);
        unc_resampling_avg = cell(1, num_grid_points);
        unc_resampling_std = cell(1, num_grid_points);
        
        % loop through grid points
        for r = 10 %1:num_rows
            for c = 10 %1:num_cols
                
                % current grid point index
                grid_point_index = (r - 1) * num_cols + c;

                % ================================================
                % extract interrogation window
                % ================================================
                % extract grid-point image co-ordinates
                X = current_result.X(r, c);
                Y = current_result.Y(r, c);
                
                % extract window size
                str = strsplit(jobfile.(['PIV' num2str(pass_number)]).winsize, ',');
                window_size_x = str2double(str{1});
                window_size_y = str2double(str{2});
                
                % extract window resolution
                str = strsplit(jobfile.(['PIV' num2str(pass_number)]).winres, ',');
                window_resolution_x = str2double(str{1});
                window_resolution_y = str2double(str{1});
                
                % identify window extents
                xmin = X - window_size_x/2 + 1;
                xmax = X + window_size_x/2;
                ymin = Y - window_size_y/2 + 1;
                ymax = Y + window_size_y/2;
                
                % find the image windows
                zone1 = im1d(max([1 ymin]):min([image_height ymax]), max([1 xmin]):min([image_width xmax]));
                zone2 = im2d(max([1 ymin]):min([image_height ymax]), max([1 xmin]):min([image_width xmax]));
                
                % minimum subtraction
                zone1 = zone1 - min(zone1(:));
                zone2 = zone2 - min(zone2(:));

                % ================================================
                % zero pad the interrogation windows
                % ================================================
                if size(zone1,1)~=window_size_y || size(zone1,2)~=window_size_x
                    w1 = zeros(window_size_y,window_size_x);
                    w1( 1+max([0 1-ymin]):window_size_y - max([0 ymax-image_height]), 1+max([0 1-xmin]):window_size_x - max([0 xmax-image_width]) ) = zone1;
                    zone1 = w1;
                end
                
                if size(zone2,1)~=window_size_y || size(zone2,2)~=window_size_x
                    w2 = zeros(window_size_y,window_size_x);
                    w2( 1+max([0 1-ymin]):window_size_y - max([0 ymax-image_height]), 1+max([0 1-xmin]):window_size_x -max([0 xmax-image_width]) ) = zone2;
                    zone2 = w2;
                end
                
                if max(zone1(:)) == 0 || max(zone2(:)) == 0
                    continue;
                end
                
                % ================================================
                %% apply windowing function to the interrogation window
                % ================================================
                
                % create window masking filter
                sfilt1 = windowmask([window_size_x window_size_y], [window_resolution_x window_resolution_y]);
                sfilt2 = windowmask([window_size_x window_size_y], [window_resolution_x window_resolution_y]);
                
                % apply the image spatial filter
                im1_sub = zone1 .* sfilt1;
                im2_sub = zone2 .* sfilt2;
                
                % ================================================
                %% identify particles
                % ================================================                
                % calculate product of images
                im_p = sqrt(im1_sub .* im2_sub); %/sqrt(max(im1_sub(:)) * max(im2_sub(:)));
                
                % set minimum intensity threshold for the image    
                intensity_threshold = prctile(im_p(:), 90);

                % identify intensity peaks using particle identification
                [p_matrix, peaks, num_p] = dynamic_threshold_segmentation_v3(im_p, intensity_threshold, 0);

                % ================================================
                %% dot sizing and centroid estimation
                % ================================================                
                % extract cross-correlation diameter
                d_c = current_result.uncertainty2D.Autod(r, c);

                % calculate particle diameter from the auto-correlation
                d_p = d_c/sqrt(2);
                
                % change size properties
                % particle sizing settings
                sizeprops_current = sizeprops;
                sizeprops_current.p_area = d_p; %0.5*d_p;
                
                % perform particle sizing
                [XYDiameter, mapsizeinfo, locxy] = particle_size_MAIN_V1(im_p, p_matrix, num_p, sizeprops_current);
                
                % extract locations of particles
                x1 = XYDiameter(:, 1);
                y1 = XYDiameter(:, 2);
                
                % extract particle diameters
                d1 = XYDiameter(:, 3);
                
                % identify nan elements
                indices = isnan(x1) | isnan(y1);
                
                % remove nan elements
                x1(indices) = [];
                y1(indices) = [];
                d1(indices) = [];
                
                % assign the positions of the peaks in the second image
                % to be the same (as the disparity is within subpixel
                % limit)
                x2 = x1;
                y2 = y1;
                
                if display_id_results
                    % ================================================
                    %% display identification results    
                    % ================================================
                    figure
                    imagesc(im_p)
                    colormap(gray)
                    caxis([0 intensity_threshold*2])
                    % maxval = 0.8 * max(im_p(:))/2;
                    % caxis([0 maxval])
                    hold on
                    plot(x1, y1, 'o')
                    
                    % plot window edge
                    x_c = size(im_p, 2)/2;
                    y_c = size(im_p, 1)/2;
                    
                    x = [x_c - window_resolution_x/2, x_c + window_resolution_x/2, x_c + window_resolution_x/2, x_c - window_resolution_x/2, x_c - window_resolution_x/2];
                    y = [y_c - window_resolution_y/2, y_c - window_resolution_y/2, y_c + window_resolution_y/2, y_c + window_resolution_y/2, y_c - window_resolution_y/2];
                    plot(x, y, '*-')
                    colorbar
                    set_axes(gca);
                    pause(0.1);
                    
                end
                
                % ================================================
                %% remove particles and perform resampling calc
                % ================================================                
                % calculate number of particles identified
                num_particles = numel(x1);
                % calculate number of particles to be removed
                num_particles_remove = round(percentage_particles_remove * num_particles);
                % calculate local variance
                tic
                [unc_resampling_avg{grid_point_index}, unc_resampling_std{grid_point_index}, unc_resampling{grid_point_index}] = calculate_variance_resampling(im1_sub, im2_sub, window_size_x/2, window_size_y/2, 0, 0, x1, y1, x2, y2, d_p*ones(size(x1)), num_particles_remove, num_resampling_trials, uncertainty_flags, 0.25*intensity_threshold, [window_size_x, window_size_y], [window_resolution_x, window_resolution_y]);                        
                toc
            end                
        end

        % =======================
        % save results to file
        % =======================
        % directory to store resampling results
        resampling_results_directory = fullfile(top_write_directory, rectype, ['Camera', num2str(camera_index), filesep], 'resampling', ['ppr' num2str(percentage_particles_remove, '%.2f') '_nrs' num2str(num_resampling_trials)]);
        mkdir_c(resampling_results_directory);
        filename = fullfile(resampling_results_directory, sprintf('unc_resampling_%05d.mat', (snapshot_index - 1) + start_frame));
        save(filename, 'unc_resampling_avg', 'unc_resampling_std', 'unc_resampling');
    end    
end

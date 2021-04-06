% This script calculates the pdf of uncertainties estimated by three
% different algorithms for experimental stereo images by jack-knife resampling.

%%
clear
close all
clc

% ============================
%% set mount directory
% ============================
if ismac
    mount_directory_a = '/Volumes/aether/';
    mount_directory_c = '/Volumes/aether_c/';
else
    mount_directory_a = '/scratch/shannon/a/aether/';
    mount_directory_c = '/scratch/shannon/c/aether';
end

% ============================
%% add paths
% ============================
restoredefaultpath;
addpath(genpath(fullfile(mount_directory_c, 'Projects/BOS/general-codes/matlab-codes/')));
setup_default_settings;
% addpath('../prana-uncertainty-average-dc-new-im-2-save-deform-cs/');
addpath('../prana/');
addpath(genpath('/scratch/shannon/c/aether/Projects/BOS/general-codes/piv-image-generation/'));
addpath('../general-codes/');
addpath('../stereo_uncertainty_codes_packaged/');

% ============================
%% experiment settings
% ============================
% seconds per frame
spf = 0.001;
% No. of frames
num_frames = 50;   
% Starting frame
fstart = 24;
% camera numbers
camera_numbers = [2, 4];
num_cameras = numel(camera_numbers);

% ====================================
%% read/write settings
% ====================================
% top level directory for this project
top_project_directory = fullfile(mount_directory_a, 'Projects/PIV_PTV_Uncertainty_Quantification/uncertainty-combination/');
% top level read directory
top_read_directory = fullfile(mount_directory_a, 'Projects/PIV_PTV_Uncertainty_Quantification/Stereo_uncertainty_work/');
% top level write directory
% top_write_directory = fullfile(top_project_directory, '/analysis/results/stereo-dataset/cam13/With_selfcal/');
top_write_directory = fullfile(top_project_directory, '/analysis/results/stereo-dataset/', ['cam' num2str(camera_numbers(1)) num2str(camera_numbers(2))], '/With_selfcal/');
mkdir_c(top_write_directory);
% Load 2d job file 
% job_settings = load(fullfile(top_write_directory, 'cam13_VR_prana_fulljob_withselfcal.mat'));
job_settings = load(fullfile(top_write_directory, ['cam' num2str(camera_numbers(1)) num2str(camera_numbers(2)) '_VR_prana_fulljob_withselfcal.mat']));
% extract calibration job file
caljobfile = job_settings.caljobfile;
% extract 2d job file
planarjob = job_settings.planarjob;

% ====================================
%% processing settings
% ====================================
% camera numbers
% camera_numbers = [1, 3];
% camera_numbers = [2, 4];
% num_cameras = numel(camera_numbers);
% stereo reconstruction type
rectype = 'Willert';
% pass number
pass_number = 4;
% array of uncertainty methods
individual_method_array = {'IM'; 'MC'; 'CS'};
num_uncertainty_methods = numel(individual_method_array);
% start frame
start_frame = 25;

% ====================================
%% monte-carlo settings
% ====================================
% num overall trials
num_trials = 1e3;

% ============================
%% resampling settings
% ============================
% resampling method names
resampling_method_name_array = {'remove-paired'; 'remove-random'; 'add-random'};
% resampling_method_name_array = {'add-random'};
num_resampling_methods = numel(resampling_method_name_array);
% number of particles to remove
percentage_particles_remove_array = 0.05:0.05:0.25; %0.025:0.025:0.25; 
num_ppr = numel(percentage_particles_remove_array);
% number of resampling trials
num_resampling_trials = 1e2;
% display_resampled_images? (true/false)
display_resampled_images = 0;

% ====================================
%% particle identification settings
% ====================================
% intensity threshold for particle identification
intensity_threshold = 1000;
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
%% load planar results
% ====================================
fprintf('Loading planar results\n');

results_planar_all = cell(1, 2);
jobfile_all = cell(1, 2);
files_im1 = cell(1, 2);
files_im2 = cell(1, 2);

% loop through cameras
for camera_index = 1:num_cameras
    fprintf('Camera: %d\n', camera_index);
    % ====================================
    %% Load data
    % ====================================
    % results directory for the current camera
    current_results_directory = fullfile(top_write_directory, rectype, ['Camera', num2str(camera_index), filesep], 'vectors');
    % load results for vectors, errors and uncertainties
    [results_planar_all{camera_index}, num_snapshots] = load_directory_data(current_results_directory, ['VR*pass' num2str(pass_number, '%d') '*.mat']);
    % extract number of grid point
    [num_rows, num_cols] = size(results_planar_all{camera_index}{1}.X);
    num_grid_points = num_rows * num_cols;
    
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
% extract job file properties
% ====================================
% extract job file
jobfile = job_settings.job1;
% extract pass settings
current_pass_settings = jobfile.(['PIV' num2str(pass_number)]);
% extract window size
[window_size_x, window_size_y] = extract_window_size(current_pass_settings);
% extract window resolution
[window_resolution_x, window_resolution_y] = extract_window_resolution(current_pass_settings);        
% extract grid resolution
[grid_resolution_x, grid_resolution_y] = extract_grid_resolution(current_pass_settings);

% ================================================
% remove edge points
% ================================================
r = 1:num_rows;
c = 1:num_cols;

% calculate number of grid points for have a window resolution
grid_point_buffer_x = 0.5 * window_resolution_x/grid_resolution_x;
grid_point_buffer_y = 0.5 * window_resolution_y/grid_resolution_y;

% update row and column list
r = r(grid_point_buffer_y:(end-grid_point_buffer_y));
c = c(grid_point_buffer_x:(end-grid_point_buffer_x));

% ====================================
%% load stereo results
% ====================================
fprintf('Loading stereo results\n');
% directory containing stereo results
stereo_results_directory = fullfile(top_write_directory, rectype, 'Camera1Camera2_3Cfields',filesep);
% load results for vectors, errors and uncertainties
[results_stereo_all, num_snapshots] = load_directory_data(fullfile(stereo_results_directory, 'vectors'), ['piv*pass_' num2str(pass_number, '%d') '*.mat']);
    
% ============================
% Extract calibration details
% ============================
% Magnification in mm/pix
scaling = job_settings.scaling.wil;
mx = scaling.xscale;
my = scaling.yscale;
% Here the results are for camera1camera3 pair
cam1 = '1';
cam2 = '3';
if strcmp(cam1,'1') && strcmp(cam2,'3')
    mz = scaling.yscale;
elseif strcmp(cam1,'2') && strcmp(cam2,'4')
    mz = scaling.xscale;
end

% ====================================
%% run resampling
% ====================================
fprintf('running monte-carlo\n');

% initialize variables to store resampling results
unc_planar_individual = cell(num_cameras, num_trials);
snr_planar_individual = cell(num_cameras, num_trials);
unc_planar_resampling = cell(num_cameras, num_trials, num_ppr, num_resampling_methods);
snr_planar_resampling = cell(num_cameras, num_trials, num_ppr, num_resampling_methods);

f_trial_all = nans(1, num_trials);
r_trial_all = nans(1, num_trials);
c_trial_all = nans(1, num_trials);

% set seed for random number generator
rng(0);
tic
for trial_index = 1:num_trials
    fprintf('trial: %d\n', trial_index);
    
    % random snapshot
    f_trial = randi(num_snapshots);

    % select a random grid point
    r_trial = randsample(r, 1);
    c_trial = randsample(c, 1);

    % loop through cameras
    for camera_index = 1:num_cameras
        % load snapshot results into memory
        current_result = results_planar_all{camera_index}{f_trial};        
        
        % ================================================
        % load deformed images
        % ================================================
        % image 1
        im1 = imread(fullfile(files_im1{camera_index}(f_trial).folder, files_im1{camera_index}(f_trial).name));
        im1d = double(flipud(im1));
        % image 2
        im2 = imread(fullfile(files_im2{camera_index}(f_trial).folder, files_im2{camera_index}(f_trial).name));
        im2d = double(flipud(im2));

        % extract image height and width
        if trial_index == 1
            [image_height, image_width] = size(im1d);
        end

        % extract grid point coordinates
        X = current_result.X(r_trial, c_trial);
        Y = current_result.Y(r_trial, c_trial);

        % extract interrogation window
        im1_sub = extract_interrogation_window(im1d, X, Y, [window_size_x, window_size_y]);
        im2_sub = extract_interrogation_window(im2d, X, Y, [window_size_x, window_size_y]);
        
        % calculate planar uncertainties
        [unc_planar_individual{camera_index, trial_index}, snr_planar_individual{camera_index, trial_index}] = calculate_planar_uncertainties(im1_sub, im2_sub, ...
                                                                    current_pass_settings, individual_method_array);

        % calculate effective particle diameter
        d_p = current_result.uncertainty2D.Autod(r_trial, c_trial)/sqrt(2);

        % ================================================
        %% loop through particle removal percentages
        % ================================================
        for particle_remove_index = 1:num_ppr
            % current percentage to remove
            percentage_particles_remove = percentage_particles_remove_array(particle_remove_index);
            % fprintf('particle remove percentage: %.2f\n', percentage_particles_remove);
            
            % ==========================
            % loop through resampling methods
            % ==========================
            for resampling_method_index = 1:num_resampling_methods
                % fprintf('resampling method index: %d\n', resampling_method_index);
                resampling_method_name = resampling_method_name_array{resampling_method_index};

                % calculate uncertainty from resampling
                [unc, snr] = calculate_uncertainty_resampling_general(im1_sub, im2_sub, d_p, sizeprops, resampling_method_name, percentage_particles_remove, ...
                                                                        num_resampling_trials, individual_method_array, 0.25*intensity_threshold, ...
                                                                        current_pass_settings, display_resampled_images);
                % copy results
                unc_planar_resampling{camera_index, trial_index, particle_remove_index, resampling_method_index} = unc;
                snr_planar_resampling{camera_index, trial_index, particle_remove_index, resampling_method_index} = snr;
            end    
        end
    end

    % ====================
    %% aggregate
    % ====================
    f_trial_all(trial_index) = f_trial;
    r_trial_all(trial_index) = r_trial;
    c_trial_all(trial_index) = c_trial;
end

% directory to store resampling results    
resampling_results_directory = fullfile(stereo_results_directory, 'resampling-new', ...
['trials=' num2str(num_trials, '%d'), '_nrs' num2str(num_resampling_trials) '-prewindow']);
mkdir_c(resampling_results_directory);
% file name to save results
filename = fullfile(resampling_results_directory, 'monte_carlo_results.mat');
% save results to file
save(filename, 'f_trial_all', 'r_trial_all', 'c_trial_all', 'unc_planar_individual', 'snr_planar_individual', ...
                'unc_planar_resampling', 'snr_planar_resampling');

toc


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
addpath('../prana-uncertainty-average-dc-new-im-2-save-deform-cs/');
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

% ====================================
%% read/write settings
% ====================================
% top level directory for this project
top_project_directory = fullfile(mount_directory_a, 'Projects/PIV_PTV_Uncertainty_Quantification/uncertainty-combination/');
% top level read directory
top_read_directory = fullfile(mount_directory_a, 'Projects/PIV_PTV_Uncertainty_Quantification/Stereo_uncertainty_work/');
% top level write directory
top_write_directory = fullfile(top_project_directory, '/analysis/results/stereo-dataset/cam13/With_selfcal/');
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
num_cameras = numel(camera_numbers);
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
%% monte-carlo settings
% ====================================
% num overall trials
num_trials = 1e3;

% ====================================
%% resampling settings
% ====================================
% number of particles to remove
% num_particles_remove_array = [5, 10, 15];
percentage_particles_remove = 0.25;
% number of resampling trials
num_resampling_trials = 1e3;
% display_resampled_images? (true/false)
display_resampled_images = 0;
% uncertainty flags for the various methods (0 = False, 1 = True)
uncertainty_flags.ppruncertainty = 0;
uncertainty_flags.miuncertainty = 1;
uncertainty_flags.mcuncertainty = 1;

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
    [num_rows, num_cols] = size(results_planar_all{camera_index}(1).X);
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
% extract window size
str = strsplit(jobfile.(['PIV' num2str(pass_number)]).winsize, ',');
window_size_x = str2double(str{1});
window_size_y = str2double(str{2});
% extract window resolution
str = strsplit(jobfile.(['PIV' num2str(pass_number)]).winres, ',');
window_resolution_x = str2double(str{1});
window_resolution_y = str2double(str{1});

% ====================================
%% load stereo results
% ====================================
fprintf('Loading stereo results\n');
% directory containing stereo results
stereo_results_directory = fullfile(top_write_directory, rectype, 'Camera1Camera2_3Cfields',filesep);
% load results for vectors, errors and uncertainties
[results_stereo_all, num_snapshots] = load_directory_data(fullfile(stereo_results_directory, 'vectors'), ['piv*pass_' num2str(pass_number, '%d') '*.mat']);
    
% load errors
fprintf('Loading all errors into memory\n');
filename = fullfile(stereo_results_directory, 'errors.mat');
errors_all = load(filename);

% load uncertainties
fprintf('Loading all uncertainties into memory\n');
filename = fullfile(stereo_results_directory, 'uncertainties.mat');
uncertainties_all = load(filename);

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
unc_planar_resampling = cell(num_cameras, num_trials);
unc_planar_resampling_avg = cell(num_cameras, num_trials);
unc_planar_resampling_std = cell(num_cameras, num_trials);

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
    r_trial = randsample(num_rows, 1);
    c_trial = randsample(num_cols, 1);

    % loop through cameras
    for camera_index = 1:num_cameras
        % load snapshot results into memory
        current_result = results_planar_all{camera_index}(f_trial);        
        
        % extract planar uncertainties
        unc_planar_individual{camera_index, trial_index} = extract_planar_uncertainties(current_result.uncertainty2D, r_trial, c_trial);
        
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

        % extract interrogation window
        im1_sub = extract_interrogation_window(current_result, im1d, [window_size_x, window_size_y], r_trial, c_trial);
        im2_sub = extract_interrogation_window(current_result, im2d, [window_size_x, window_size_y], r_trial, c_trial);

        % calculate effective particle diameter
        d_p = current_result.uncertainty2D.Autod(r_trial, c_trial)/sqrt(2);

        % identify particles
        [x1, y1, x2, y2] = identify_particles(im1_sub, im2_sub, d_p, sizeprops, display_id_results);
        
        % ================================================
        %% calculate uncertainty by resampling particles
        % ================================================                
        % calculate number of particles identified
        num_particles = numel(x1);

        % calculate number of particles to be removed
        num_particles_remove = round(percentage_particles_remove * num_particles);
        
        % calculate uncertainty
        [unc_planar_resampling_avg{camera_index, trial_index}, unc_planar_resampling_std{camera_index, trial_index}, unc_planar_resampling{camera_index, trial_index}] = calculate_uncertainty_resampling(im1_sub, im2_sub, x1, y1, x2, y2, ...
                                d_p*ones(size(x1)), num_particles_remove, num_resampling_trials, uncertainty_flags, mean(im1_sub(:)), ... %0.25*intensity_threshold, ...
                                [window_size_x, window_size_y], [window_resolution_x, window_resolution_y], display_resampled_images);
    end

    % ====================
    %% aggregate
    % ====================
    f_trial_all(trial_index) = f_trial;
    r_trial_all(trial_index) = r_trial;
    c_trial_all(trial_index) = c_trial;
end

% directory to store resampling results    
resampling_results_directory = fullfile(stereo_results_directory, 'resampling', ...
['trials=' num2str(num_trials, '%d')], ['ppr' num2str(percentage_particles_remove, '%.2f') '_nrs' num2str(num_resampling_trials) '-prewindow']);
mkdir_c(resampling_results_directory);
% file name to save results
filename = fullfile(resampling_results_directory, 'monte_carlo_results.mat');
% save results to file
save(filename, 'f_trial_all', 'r_trial_all', 'c_trial_all', 'unc_planar_individual', ...
                'unc_planar_resampling_avg', 'unc_planar_resampling_std', 'unc_planar_resampling');
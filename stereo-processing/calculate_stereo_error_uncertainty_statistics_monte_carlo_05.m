% This script calculates the stereo uncertainty of individual and combined
% schemes, by propagation through the stereo measurement chain.

%%
clear
close all
clc

restoredefaultpath;
addpath(genpath('/scratch/shannon/c/aether/Projects/BOS/general-codes/matlab-codes/'));
addpath('../prana-uncertainty-average-dc-new-im-2-save-deform-cs/');
addpath('../general-codes/');
addpath('../stereo_uncertainty_codes_packaged/');
addpath('../histogram_distance/') 

setup_default_settings;
% dbstop if error

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
% top level read directory
top_read_directory = '/scratch/shannon/a/aether/Projects/PIV_PTV_Uncertainty_Quantification/Stereo_uncertainty_work/';
% top level write directory
% top_write_directory = '/scratch/shannon/a/aether/Projects/PIV_PTV_Uncertainty_Quantification/uncertainty-combination/analysis/results/stereo-dataset/cam13/With_selfcal/';
top_write_directory = fullfile('/scratch/shannon/a/aether/Projects/PIV_PTV_Uncertainty_Quantification/uncertainty-combination/', '/analysis/results/stereo-dataset/', ['cam' num2str(camera_numbers(1)) num2str(camera_numbers(2))], '/With_selfcal/');

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
camera_numbers = [1, 3];
num_cameras = numel(camera_numbers);
% stereo reconstruction type
rectype = 'Willert';
% pass number
pass_number = 4;
% array of uncertainty methods
individual_method_array = {'IM'; 'MC'; 'CS'};
num_individual_methods = numel(individual_method_array);
% name of the uncertainty methods in prana
uncertainty_method_prana_names = {'Uim'; 'MC'; 'Ucs'};
% start frame
start_frame = 25;
% component names
component_names = {'x'; 'y'; 'z'};
% number of velocity component_names
num_components = 3;

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
% num overall trials
num_trials = 1e3;

% ========================
%% weight settings
% ========================
% number of bins for weights
num_bins_w = 30;
% bins for weights
bins_w = linspace(0, 1, num_bins_w);
% metric type ('rms', 'iqr', 'std')
metric_name = 'iqr';

% ===============================
%% statistical analysis settings
% ===============================
% minimum allowable error (pix.)
error_min = 1e-4;
% maximum allowable error (pix.)
error_max = [0.5, 0.5, 1.5, 1.5];

% number of bins for histogram
num_bins = 30;
% bins for histograms
bins = cell(1, num_components+1);
for component_index = 1:(num_components+1)
    bins{component_index} = linspace(error_min, error_max(component_index), num_bins);
end

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
end

% ====================================
% extract job file properties
% ====================================
% extract job file
jobfile = job_settings.job1;

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

% ===============================
%% load resampling results
% ===============================
fprintf('loading resampling results\n');
% directory to load results for this case
current_read_directory = fullfile(stereo_results_directory, 'resampling-new', ...
['trials=' num2str(num_trials, '%d') '_nrs' num2str(num_resampling_trials) '-prewindow']);

% load monte carlo results
load(fullfile(current_read_directory, 'monte_carlo_results.mat'));

% ===============================
%% load weights and combined uncertainty
% ===============================
fprintf('loading weights and combined uncertainty\n');
load(fullfile(current_read_directory, 'combined_uncertainties_3c.mat'));

% ===============================
%% initialize variables 
% ===============================
fprintf('aggregating results\n');
% err_all = nans(3, num_trials);
err_all = struct;
unc_stereo_individual_all = struct;
unc_stereo_combined_all = struct;
weights_all = cell(1, num_resampling_methods);

for component_index = 1:num_components
    component_name = component_names{component_index};
    err_all.(component_name) = nans(num_trials, 1);

    unc_stereo_individual_all.(component_name) = nans(num_trials, num_individual_methods);
    unc_stereo_combined_all.(component_name) = nans(num_trials, num_resampling_methods);

    for method_index = 1:num_resampling_methods
        weights_all{method_index}.(component_name) = nans(num_trials, num_individual_methods);
    end
end

% unc_stereo_individual_all = {nans(3, num_trials); nans(3, num_trials); nans(3, num_trials)};
% unc_stereo_resampled_all = {nans(3, num_trials * num_resampling_trials); nans(3, num_trials * num_resampling_trials); nans(3, num_trials * num_resampling_trials)};
% unc_stereo_combined_all = {nans(3, num_trials); nans(3, num_trials); nans(3, num_trials)};
% weights_all = cell(num_resampling_methods, num_individual_methods);

% unc_planar_individual_all = {nans(3, num_trials), nans(3, num_trials); nans(3, num_trials), nans(3, num_trials)};
% unc_planar_resampled_all = {nans(3, num_trials * num_resampling_trials), nans(3, num_trials * num_resampling_trials); ...
%                                 nans(3, num_trials * num_resampling_trials), nans(3, num_trials * num_resampling_trials)};


unc_planar_individual_all = {nans(2, num_trials), nans(2, num_trials), nans(2, num_trials); nans(2, num_trials), nans(2, num_trials), nans(2, num_trials)};
unc_planar_resampled_all = {nans(2, num_trials * num_resampling_trials), nans(2, num_trials * num_resampling_trials), nans(2, num_trials * num_resampling_trials); ...
                            nans(2, num_trials * num_resampling_trials), nans(2, num_trials * num_resampling_trials), nans(2, num_trials * num_resampling_trials)};

% ===============================
% loop through trials
% ===============================
for trial_index = 1:num_trials
    fprintf('trial: %d\n', trial_index);
    % extract snapshot for this trial
    f_trial = f_trial_all(trial_index);

    % extract grid points for this trial
    r_trial = r_trial_all(trial_index);
    c_trial = c_trial_all(trial_index);
        
    % aggregate planar uncertainties - individual
    for camera_index = 1:num_cameras
        for method_index = 1:num_individual_methods
            for component_index = 1:2
                unc_planar_individual_all{camera_index, method_index}(component_index, trial_index) = real(unc_planar_individual{camera_index, trial_index}{method_index}.(component_names{component_index}));
            end
        end
    end

    % % aggregate planar uncertainties - individual, resampled
    % start_index = (trial_index - 1) * num_resampling_trials + 1;
    % stop_index = trial_index * num_resampling_trials;
    % for camera_index = 1:num_cameras
    %     for method_index = 1:num_individual_methods
    %         for resampling_method_index = 1:num_resampling_methods
    %             for particle_remove_index = 1:num_ppr
    %                 for component_index = 1:2
    %                     unc_planar_resampled_all{camera_index, particle_remove_index, resampling_method_index, method_index}(component_index, start_index:stop_index) = real(unc_planar_resampling{camera_index, trial_index, particle_remove_index, resampling_method_index}.([lower(individual_method_array{method_index}) component_names{component_index}]));
    %                 end
    %             end
    %         end
    %     end
    % end
    
    % aggregate stereo errors
    err_all.x(trial_index) = errors_all.err_u(r_trial, c_trial, f_trial);
    err_all.y(trial_index) = errors_all.err_v(r_trial, c_trial, f_trial);
    err_all.z(trial_index) = errors_all.err_w(r_trial, c_trial, f_trial);    
    
    for component_index = 1:num_components
        component_name = component_names{component_index};
        
        % aggregate stereo uncertainties - individual
        for method_index = 1:num_individual_methods
            unc_stereo_individual_all.(component_name)(trial_index, method_index) = real(unc_stereo_individual{trial_index}{method_index}.(component_name));
        end

        % aggregate stereo uncertainties - combined
        for method_index = 1:num_resampling_methods
            unc_stereo_combined_all.(component_name)(trial_index, method_index) = real(unc_stereo_combined{trial_index}{method_index}.(component_name));
        end

        % aggregate weights    
        for combination_method_index = 1:num_resampling_methods
            for uncertainty_method_index = 1:num_individual_methods
                weights_all{combination_method_index}.(component_name)(trial_index, uncertainty_method_index) = weights{trial_index}{combination_method_index}.(component_names{component_index})(uncertainty_method_index);
            end
        end    
    end
end


% extract valid measurements
[valid_trials, err_all_valid, unc_stereo_individual_all_valid, unc_stereo_combined_all_valid] = extract_valid_errors_uncertainties(err_all, unc_stereo_individual_all, unc_stereo_combined_all, ...
                                                                                        component_names, error_min*ones(1, num_components), error_max(1:num_components));
num_valid_trials = sum(double(valid_trials));

% ===============================
%% calculate rms of uncertainties
% ===============================
fprintf('calculating rms\n');
err_rms = nans(1, num_components+1);
unc_indiv_rms = nans(num_components, num_individual_methods);
unc_comb_rms = nans(num_components, num_resampling_methods);
for component_index = 1:num_components
    component_name = component_names{component_index};
    % errors
    err_rms(component_index) = rms(err_all_valid.(component_name));
    % uncertainties - individual
    unc_indiv_rms(component_index, :) = rms(unc_stereo_individual_all_valid.(component_name), 1);
    % uncertainties - combined
    unc_comb_rms(component_index, :) = rms(unc_stereo_combined_all_valid.(component_name), 1);
end

% ===============================
%% calculate pdf of weights
% ===============================
pdf_weights = calculate_pdf_weights_stereo(weights_all, bins_w, valid_trials, num_individual_methods, num_resampling_methods, component_names);

% ===================================
%% estimate error from uncertainties
% ===================================
% seed random number generator
rng(100);
% -----------------------
% individual methods
% -----------------------
num_methods = num_individual_methods;
% err_est_valid_individual = cell(num_methods, num_components+1);

for component_index = 1:num_components
    component_name = component_names{component_index};
    for method_index = 1:num_methods
        % estimate error from uncertainty
        err_est_current = randn(num_valid_trials, 1) .* unc_stereo_individual_all_valid.(component_name)(:, method_index); 
        err_est_valid_individual.(component_name)(:, method_index) = err_est_current;
    end
end

% -----------------------
% combined methods
% -----------------------
num_methods = num_resampling_methods;
% err_est_valid_combined = cell(num_methods, num_components+1);

for component_index = 1:num_components
    component_name = component_names{component_index};
    for method_index = 1:num_methods
        % estimate error from uncertainty
        err_est_current = randn(num_valid_trials, 1) .* unc_stereo_combined_all_valid.(component_name)(:, method_index); 
        err_est_valid_combined.(component_name)(:, method_index) = err_est_current;
    end
end

% ==================
% save results
% ==================
fprintf('saving results to file\n');
filename = fullfile(current_read_directory, 'stereo_errors_uncertainties_statistics_05.mat');
save(filename, 'err_all', 'unc_planar_individual_all', 'unc_stereo_individual_all', 'unc_stereo_combined_all', 'valid_trials', ...
                'err_all_valid', 'unc_stereo_individual_all_valid', 'unc_stereo_combined_all_valid', ...                
                'err_rms', 'unc_indiv_rms', 'unc_comb_rms', 'weights_all', 'pdf_weights', ...
                'err_est_valid_individual', 'err_est_valid_combined');



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
% bins for histogram of weights
num_bins_weights = 25;
bins_weights = linspace(0, 1, num_bins_weights);
% small value to replace zero
small_value = 1e-10;

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
%% aggregate results
% ===============================
fprintf('aggregating results\n');
err_all = nans(3, num_trials);
unc_stereo_individual_all = {nans(3, num_trials); nans(3, num_trials); nans(3, num_trials)};
unc_stereo_resampled_all = {nans(3, num_trials * num_resampling_trials); nans(3, num_trials * num_resampling_trials); nans(3, num_trials * num_resampling_trials)};
unc_stereo_combined_all = {nans(3, num_trials); nans(3, num_trials); nans(3, num_trials)};
weights_all = cell(num_resampling_methods, num_individual_methods);

% unc_planar_individual_all = {nans(3, num_trials), nans(3, num_trials); nans(3, num_trials), nans(3, num_trials)};
% unc_planar_resampled_all = {nans(3, num_trials * num_resampling_trials), nans(3, num_trials * num_resampling_trials); ...
%                                 nans(3, num_trials * num_resampling_trials), nans(3, num_trials * num_resampling_trials)};


unc_planar_individual_all = {nans(2, num_trials), nans(2, num_trials), nans(2, num_trials); nans(2, num_trials), nans(2, num_trials), nans(2, num_trials)};
unc_planar_resampled_all = {nans(2, num_trials * num_resampling_trials), nans(2, num_trials * num_resampling_trials), nans(2, num_trials * num_resampling_trials); ...
                            nans(2, num_trials * num_resampling_trials), nans(2, num_trials * num_resampling_trials), nans(2, num_trials * num_resampling_trials)};

% loop through trials
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
    % err_all(1, trial_index) = abs(errors_all.err_u(r_trial, c_trial, f_trial));
    % err_all(2, trial_index) = abs(errors_all.err_v(r_trial, c_trial, f_trial));
    % err_all(3, trial_index) = abs(errors_all.err_w(r_trial, c_trial, f_trial));
    err_all(1, trial_index) = errors_all.err_u(r_trial, c_trial, f_trial);
    err_all(2, trial_index) = errors_all.err_v(r_trial, c_trial, f_trial);
    err_all(3, trial_index) = errors_all.err_w(r_trial, c_trial, f_trial);
     
    % aggregate stereo uncertainties - individual
    for method_index = 1:num_individual_methods
        for component_index = 1:num_components
            unc_stereo_individual_all{method_index}(component_index, trial_index) = real(unc_stereo_individual{trial_index}{method_index}.(component_names{component_index}));
        end
    end

    % % aggregate stereo uncertainties - individual, resampled
    % start_index = (trial_index - 1) * num_resampling_trials + 1;
    % stop_index = trial_index * num_resampling_trials;
    % for method_index = 1:num_individual_methods
    %     unc_stereo_resampled_all{method_index}(1, start_index:stop_index) = real(unc_stereo_resampled{trial_index}{method_index}.x);
    %     unc_stereo_resampled_all{method_index}(2, start_index:stop_index) = real(unc_stereo_resampled{trial_index}{method_index}.y);
    %     unc_stereo_resampled_all{method_index}(3, start_index:stop_index) = real(unc_stereo_resampled{trial_index}{method_index}.z);
    % end

    % aggregate stereo uncertainties - combined
    for method_index = 1:num_resampling_methods
        for component_index = 1:num_components
            unc_stereo_combined_all{method_index}(component_index, trial_index) = real(unc_stereo_combined{trial_index}{method_index}.(component_names{component_index}));
        end
    end

    % aggregate weights    
    for combination_method_index = 1:num_resampling_methods
        for uncertainty_method_index = 1:num_individual_methods
            for component_index = 1:num_components
                weights_all{combination_method_index, uncertainty_method_index}(component_index, trial_index) = weights{trial_index}{combination_method_index}.(component_names{component_index})(uncertainty_method_index);
            end
        end
    end
end

% ===============================
%% remove invalid measurements
% ===============================
fprintf('removing invalid measurements\n');

% errors
err_all_valid = cell(1, num_components+1);
err_all_valid{num_components+1} = [];
err_all = real(err_all);
for component_index = 1:num_components
    % err_all_valid{component_index} = remove_invalid_measurements(abs(err_all(component_index, :)), error_min, error_max(component_index));
    err_all_valid{component_index} = remove_invalid_measurements(err_all(component_index, :), -error_max(component_index), error_max(component_index));
    err_all_valid{num_components+1} = [err_all_valid{num_components+1}, err_all_valid{component_index}];
end

% planar uncertainties - individual
% unc_planar_individual_all_valid = cell(num_individual_methods, 2, num_cameras);
unc_planar_individual_all_valid = cell(num_cameras, num_individual_methods);
% valid_indices_individual = cell(num_resampling_methods, num_components);
for camera_index = 1:num_cameras
    for method_index = 1:num_individual_methods
        unc_planar_individual_all_valid{camera_index, method_index} = cell(1, 2);
        for component_index = 1:2
            % unc_planar_individual_all_valid{method_index, component_index, camera_index} = remove_invalid_measurements(unc_planar_individual_all{camera_index, method_index}(component_index, :), error_min, error_max(component_index));            
            unc_planar_individual_all_valid{camera_index, method_index}{component_index} = remove_invalid_measurements(unc_planar_individual_all{camera_index, method_index}(component_index, :), error_min, error_max(component_index));            
        end    
    end
end

% % planar uncertainties - individual, resampled
% % unc_planar_resampled_all_valid = cell(num_individual_methods, 2, num_cameras);
% unc_planar_resampled_all_valid = cell(num_cameras, num_individual_methods);
% % valid_indices_individual = cell(num_resampling_methods, num_components);
% for camera_index = 1:num_cameras
%     for method_index = 1:num_individual_methods
%         unc_planar_resampled_all_valid{camera_index, method_index} = cell(1, 2);
%         for component_index = 1:2
%             % unc_planar_resampled_all_valid{method_index, component_index, camera_index} = remove_invalid_measurements(unc_planar_resampled_all{camera_index, method_index}(component_index, :), error_min, error_max(component_index));            
%             unc_planar_resampled_all_valid{camera_index, method_index}{component_index} = remove_invalid_measurements(unc_planar_resampled_all{camera_index, method_index}(component_index, :), error_min, error_max(component_index));            
%         end    
%     end
% end

% stereo uncertainties - individual
unc_stereo_individual_all_valid = cell(num_individual_methods, num_components+1);
valid_indices_individual = cell(num_resampling_methods, num_components);
for method_index = 1:num_individual_methods
    unc_stereo_individual_all_valid{method_index, num_components+1} = [];
    for component_index = 1:num_components
        [unc_stereo_individual_all_valid{method_index, component_index}, valid_indices_individual{method_index, component_index}] = remove_invalid_measurements(unc_stereo_individual_all{method_index}(component_index, :), error_min, error_max(component_index));
        unc_stereo_individual_all_valid{method_index, num_components+1} = [unc_stereo_individual_all_valid{method_index, num_components+1}, ...
                                                                        unc_stereo_individual_all_valid{method_index, component_index}];
    end    
end

% % stereo uncertainties - individual, resampled
% unc_stereo_resampled_all_valid = cell(num_individual_methods, num_components+1);
% for method_index = 1:num_individual_methods
%     unc_stereo_resampled_all_valid{method_index, num_components+1} = [];
%     % loop through co-ordinates
%     for component_index = 1:num_components
%         % identify valid indices
%         unc_stereo_resampled_all_valid{method_index, component_index} = remove_invalid_measurements(unc_stereo_resampled_all{method_index}(component_index, :), error_min, error_max(component_index));
%         unc_stereo_resampled_all_valid{method_index, num_components+1} = [unc_stereo_resampled_all_valid{method_index, num_components+1}, ...
%                                                                         unc_stereo_resampled_all_valid{method_index, component_index}];
%     end    
% end

% stereo uncertainties - combined
unc_stereo_combined_all_valid = cell(num_resampling_methods, num_components+1);
valid_indices_combined = cell(num_resampling_methods, num_components);
for method_index = 1:num_resampling_methods
    unc_stereo_combined_all_valid{method_index, num_components+1} = [];
    for component_index = 1:num_components
        [unc_stereo_combined_all_valid{method_index, component_index}, valid_indices_combined{method_index, component_index}] = remove_invalid_measurements(unc_stereo_combined_all{method_index}(component_index, :), error_min, error_max(component_index));
        unc_stereo_combined_all_valid{method_index, num_components+1} = [unc_stereo_combined_all_valid{method_index, num_components+1}, ...
                                                                        unc_stereo_combined_all_valid{method_index, component_index}];
    end
end

% ===============================
%% calculate pdf and cdf of uncertainties
% ===============================
fprintf('calculating pdf and cdf\n');

% errors
pdf_err = cell(1, num_components+1);
cdf_err = cell(1, num_components+1);
for component_index = 1:(num_components+1)
    pdf_err{component_index} = histcounts(abs(err_all_valid{component_index}), bins{component_index}, 'normalization', 'pdf');
    cdf_err{component_index} = histcounts(abs(err_all_valid{component_index}), bins{component_index}, 'normalization', 'cdf');
end

% uncertainties - individual
pdf_unc_individual = cell(num_individual_methods, num_components+1);
cdf_unc_individual = cell(num_individual_methods, num_components+1);
for method_index = 1:num_individual_methods
    for component_index = 1:(num_components+1)
        pdf_unc_individual{method_index, component_index} = histcounts(unc_stereo_individual_all_valid{method_index, component_index}, bins{component_index}, 'normalization', 'pdf');
        cdf_unc_individual{method_index, component_index} = histcounts(unc_stereo_individual_all_valid{method_index, component_index}, bins{component_index}, 'normalization', 'cdf');
    end
end

% % uncertainties - individual, resampled
% pdf_unc_resampled = cell(num_individual_methods, num_components+1);
% cdf_unc_resampled = cell(num_individual_methods, num_components+1);
% for method_index = 1:num_individual_methods
%     for component_index = 1:(num_components+1)
%         pdf_unc_resampled{method_index, component_index} = histcounts(unc_stereo_resampled_all_valid{method_index, component_index}, bins{component_index}, 'normalization', 'pdf');
%         cdf_unc_resampled{method_index, component_index} = histcounts(unc_stereo_resampled_all_valid{method_index, component_index}, bins{component_index}, 'normalization', 'cdf');
%     end
% end

% uncertainties - combined
pdf_unc_combined = cell(num_resampling_methods, num_components+1);
cdf_unc_combined = cell(num_resampling_methods, num_components+1);
for method_index = 1:num_resampling_methods
    for component_index = 1:(num_components+1)
        pdf_unc_combined{method_index, component_index} = histcounts(unc_stereo_combined_all_valid{method_index, component_index}, bins{component_index}, 'normalization', 'pdf');
        cdf_unc_combined{method_index, component_index} = histcounts(unc_stereo_combined_all_valid{method_index, component_index}, bins{component_index}, 'normalization', 'cdf');
    end
end

% ===============================
%% calculate rms of uncertainties
% ===============================
fprintf('calculating rms\n');

% errors
rms_err = nans(1, num_components+1);
for component_index = 1:(num_components+1)
    rms_err(component_index) = rms(err_all_valid{component_index});
end

% uncertainties - individual
rms_unc_individual = nans(num_individual_methods, num_components+1);
for method_index = 1:num_individual_methods
    for component_index = 1:(num_components+1)
        rms_unc_individual(method_index, component_index) = rms(unc_stereo_individual_all_valid{method_index, component_index});
    end
end

% % uncertainties - individual, resampled
% rms_unc_resampled = nans(num_individual_methods, num_components+1);
% for method_index = 1:num_individual_methods
%     for component_index = 1:(num_components+1)
%         rms_unc_resampled(method_index, component_index) = rms(unc_stereo_resampled_all_valid{method_index, component_index});
%     end
% end

% uncertainties - combined
rms_unc_combined = nans(num_resampling_methods, num_components+1);
for method_index = 1:num_resampling_methods
    for component_index = 1:(num_components+1)
        rms_unc_combined(method_index, component_index) = rms(unc_stereo_combined_all_valid{method_index, component_index});
    end
end

% ===============================
%% calculate pdf of weights
% ===============================
fprintf('calculating pdf of weights\n');

pdf_weights = cell(num_resampling_methods, num_individual_methods);
weights_all_valid = cell(num_resampling_methods, num_individual_methods);
% loop through combination methods
for combination_method_index = 1:num_resampling_methods
    % loop through individual methods    
    for uncertainty_method_index = 1:num_individual_methods
        % allocate memory
        pdf_weights{combination_method_index, uncertainty_method_index} = nans(3, num_bins_weights - 1);
        weights_all_valid{combination_method_index, uncertainty_method_index} = cell(1, num_components);
        % loop through u, v, and w
        for component_index = 1:num_components
            % extract weights for valid measurements
            weights_all_valid{combination_method_index, uncertainty_method_index}{component_index} = weights_all{combination_method_index, uncertainty_method_index}(component_index, valid_indices_combined{combination_method_index, component_index});
            % calculate pdf of weights 
            [pdf_weights{combination_method_index, uncertainty_method_index}(component_index, :), ~] = histcounts(weights_all_valid{combination_method_index, uncertainty_method_index}{component_index}, bins_weights, 'normalization', 'pdf');
        end
    end
end

% ===================================
%% estimate error from uncertainties
% ===================================
% seed random number generator
rng(100);

% -----------------------
% individual methods
% -----------------------
num_methods = num_individual_methods;

err_est_valid_individual = cell(num_methods, num_components+1);

for component_index = 1:(num_components+1)
    for method_index = 1:num_methods
        % calculate number of valid elements
        num_valid_trials = numel(unc_stereo_individual_all_valid{method_index, component_index});
        % estimate error from uncertainty
        err_est_current = randn(1, num_valid_trials) .* unc_stereo_individual_all_valid{method_index, component_index};
        % remove invalid measurements
        % err_est_valid_individual{method_index, component_index} = remove_invalid_measurements(abs(err_est_current), error_min, error_max(component_index));
        err_est_valid_individual{method_index, component_index} = remove_invalid_measurements(err_est_current, -error_max(component_index), error_max(component_index));
    end
end

% -----------------------
% combined methods
% -----------------------
num_methods = num_resampling_methods;

err_est_valid_combined = cell(num_methods, num_components+1);

for component_index = 1:(num_components+1)
    for method_index = 1:num_methods
        % calculate number of valid elements
        num_valid_trials = numel(unc_stereo_combined_all_valid{method_index, component_index});
        % estimate error from uncertainty
        err_est_current = randn(1, num_valid_trials) .* unc_stereo_combined_all_valid{method_index, component_index};
        % remove invalid measurements
        % err_est_valid_combined{method_index, component_index} = remove_invalid_measurements(abs(err_est_current), error_min, error_max(component_index));
        err_est_valid_combined{method_index, component_index} = remove_invalid_measurements(err_est_current, -error_max(component_index), error_max(component_index));
    end
end

% ===================================
%% calculate statistics for error estimated from uncertainties
% ===================================

% -----------------------
% error
% -----------------------
h_err = cell(num_components+1);

for component_index = 1:(num_components+1)
    % calculate histogram of estimated error
    h_err{component_index} = histcounts(abs(err_all_valid{component_index}), bins{component_index});
    % replace zero values with small finite number
    indices = h_err{component_index} == 0;
    h_err{component_index}(indices) = small_value;    
end

% -----------------------
% individual methods
% -----------------------
num_methods = num_individual_methods;
% allocate memory
pdf_err_est_individual = cell(num_components+1);
cdf_err_est_individual = cell(num_components+1);
h_err_est_individual = cell(num_components+1);
d_err_est_individual = cell(num_components+1);

% loop through component_names
for component_index = 1:(num_components+1)
    % initialize arrays
    pdf_err_est_individual{component_index} = nans(num_methods, numel(bins{component_index})-1);
    cdf_err_est_individual{component_index} = nans(num_methods, numel(bins{component_index})-1);
    h_err_est_individual{component_index} = nans(num_methods, numel(bins{component_index})-1);
    d_err_est_individual{component_index} = nans(num_methods, 1);

    % loop through methods
    for method_index = 1:num_methods
        % calculate pdf of estimated error
        pdf_err_est_individual{component_index}(method_index, :) = histcounts(abs(err_est_valid_individual{method_index, component_index}), bins{component_index}, 'normalization', 'pdf');
        % calculate cdf of estimated error
        cdf_err_est_individual{component_index}(method_index, :) = histcounts(abs(err_est_valid_individual{method_index, component_index}), bins{component_index}, 'normalization', 'cdf');
        % calculate histogram of estimated error
        h_err_est_individual{component_index}(method_index, :) = histcounts(abs(err_est_valid_individual{method_index, component_index}), bins{component_index});
        % replace zero values with small finite number
        indices = h_err_est_individual{component_index}(method_index, :) == 0;
        h_err_est_individual{component_index}(method_index, indices) = small_value;
        % calculate histogram distance measure
        d_err_est_individual{component_index}(method_index) = pdist2(h_err{component_index}, h_err_est_individual{component_index}(method_index, :), @total_variation_distance);
    end
end

% -----------------------
% combined methods
% -----------------------
num_methods = num_resampling_methods;

% allocate memory
pdf_err_est_combined = cell(num_components+1);
cdf_err_est_combined = cell(num_components+1);
h_err_est_combined = cell(num_components+1);
d_err_est_combined = cell(num_components+1);

for component_index = 1:(num_components+1)
    % initialize arrays
    pdf_err_est_combined{component_index} = nans(num_methods, numel(bins{component_index})-1);
    cdf_err_est_combined{component_index} = nans(num_methods, numel(bins{component_index})-1);
    h_err_est_combined{component_index} = nans(num_methods, numel(bins{component_index})-1);
    d_err_est_combined{component_index} = nans(num_methods, 1);

    % loop through methods
    for method_index = 1:num_methods
        % calculate pdf of estimated error
        pdf_err_est_combined{component_index}(method_index, :) = histcounts(abs(err_est_valid_combined{method_index, component_index}), bins{component_index}, 'normalization', 'pdf');
        % calculate cdf of estimated error
        cdf_err_est_combined{component_index}(method_index, :) = histcounts(abs(err_est_valid_combined{method_index, component_index}), bins{component_index}, 'normalization', 'cdf');
        % calculate histogram of estimated error
        h_err_est_combined{component_index}(method_index, :) = histcounts(abs(err_est_valid_combined{method_index, component_index}), bins{component_index});
        % replace zero values with small finite number
        indices = h_err_est_combined{component_index}(method_index, :) == 0;
        h_err_est_combined{component_index}(method_index, indices) = small_value;
        % calculate histogram distance measure
        d_err_est_combined{component_index}(method_index) = pdist2(h_err{component_index}, h_err_est_combined{component_index}(method_index, :), @total_variation_distance);
    end    
end

% ==================
% save results
% ==================
fprintf('saving results to file\n');
filename = fullfile(current_read_directory, 'stereo_errors_uncertainties_statistics_03.mat');
% save(filename, 'err_all', 'unc_planar_individual_all', 'unc_planar_resampled_all', ...
%                 'unc_stereo_individual_all', 'unc_stereo_resampled_all', 'unc_stereo_combined_all', 'valid_indices_combined', ...
%                 'err_all_valid', 'unc_planar_individual_all_valid', 'unc_planar_resampled_all_valid', ...
%                 'unc_stereo_individual_all_valid', 'unc_stereo_resampled_all_valid', 'unc_stereo_combined_all_valid', ...
%                 'pdf_err', 'cdf_err', 'pdf_unc_individual', 'cdf_unc_individual', 'pdf_unc_resampled', 'cdf_unc_resampled', ...
%                 'pdf_unc_combined', 'cdf_unc_combined', ...
%                 'rms_err', 'rms_unc_individual', 'rms_unc_resampled', 'rms_unc_combined', ...
%                 'weights_all', 'weights_all_valid', 'pdf_weights', ...
%                 'err_est_valid_individual', 'pdf_err_est_individual', 'cdf_err_est_individual', 'h_err_est_individual', 'd_err_est_individual', ...
%                 'err_est_valid_combined', 'pdf_err_est_combined', 'cdf_err_est_combined', 'h_err_est_combined', 'd_err_est_combined', 'h_err');
save(filename, 'err_all', 'unc_planar_individual_all', ...
                'unc_stereo_individual_all', 'unc_stereo_combined_all', 'valid_indices_combined', ...
                'err_all_valid', 'unc_planar_individual_all_valid', ...
                'unc_stereo_individual_all_valid', 'unc_stereo_combined_all_valid', ...
                'pdf_err', 'cdf_err', 'pdf_unc_individual', 'cdf_unc_individual', ...
                'pdf_unc_combined', 'cdf_unc_combined', ...
                'rms_err', 'rms_unc_individual', 'rms_unc_combined', ...
                'weights_all', 'weights_all_valid', 'pdf_weights', ...
                'err_est_valid_individual', 'pdf_err_est_individual', 'cdf_err_est_individual', 'h_err_est_individual', 'd_err_est_individual', ...
                'err_est_valid_combined', 'pdf_err_est_combined', 'cdf_err_est_combined', 'h_err_est_combined', 'd_err_est_combined', 'h_err');



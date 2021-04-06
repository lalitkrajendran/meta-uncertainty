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
uncertainty_method_array = {'IM'; 'MC'; 'CS'};
num_uncertainty_methods = numel(uncertainty_method_array);
% name of the uncertainty methods in prana
uncertainty_method_prana_names = {'Uim'; 'MC'; 'Ucs'};
% start frame
start_frame = 25;

% ====================================
%% resampling settings
% ====================================
% number of particles to remove
% num_particles_remove_array = [5, 10, 15];
percentage_particles_remove = 0.25;
% number of resampling trials
num_resampling_trials = 1e3;
% num overall trials
num_trials = 1e3;
% uncertainty flags for the various methods (0 = False, 1 = True)
uncertainty_flags.ppruncertainty = 0;
uncertainty_flags.miuncertainty = 0;
uncertainty_flags.mcuncertainty = 1;

% ====================================
% combination method settings
% ====================================
% array of uncertainty combination model to be used
% 'unweighted', 'global-weight-var', 'global-weight-std',
% 'local-weight-var', 'local-weight-std'
combination_methods = {'unwt'; 'var'; 'entropy'}; %'prob'};
num_combination_methods = numel(combination_methods);

% ===============================
%% analysis settings
% ===============================
% minimum allowable error (pix.)
error_min = 1e-4;
% maximum allowable error (pix.)
error_max = [0.5, 0.5, 1.5];

% number of bins for histogram
num_bins = 30;
% bins for histograms
bins = cell(1, 4);
for i = 1:3
    bins{i} = linspace(error_min, error_max(i), num_bins);
end
bins{4} = bins{3};

% number of bins for the coarse histogram
num_bins_coarse = 8;
% bins for histogram of weights
bins_weights = linspace(0, 1, 25);

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
current_read_directory = fullfile(stereo_results_directory, 'resampling', ...
['trials=' num2str(num_trials, '%d')], ['ppr' num2str(percentage_particles_remove, '%.2f') '_nrs' num2str(num_resampling_trials)]);

% load monte carlo results
load(fullfile(current_read_directory, 'monte_carlo_results.mat'));

% ===============================
%% load weights and combined uncertainty
% ===============================
fprintf('loading weights and combined uncertainty\n');
load(fullfile(current_read_directory, 'combined_uncertainties.mat'));

% ===============================
%% load stereo errors and uncertainties
% ===============================
fprintf('loading stereo errors and uncertainties\n');
load(fullfile(current_read_directory, 'stereo_errors_uncertainties.mat'));

% ===============================
%% aggregate results
% ===============================
fprintf('aggregating results\n');
err_all = nans(3, num_trials);
unc_stereo_individual_all = {nans(3, num_trials); nans(3, num_trials); nans(3, num_trials);};
unc_stereo_combined_all = {nans(3, num_trials); nans(3, num_trials); nans(3, num_trials);};

for trial_index = 1:num_trials
    fprintf('trial: %d\n', trial_index);

    % errors
    err_all(1, trial_index) = abs(errors{trial_index}.x);
    err_all(2, trial_index) = abs(errors{trial_index}.y);
    err_all(3, trial_index) = abs(errors{trial_index}.z);
     
    % uncertainties - individual
    for method_index = 1:num_uncertainty_methods
        unc_stereo_individual_all{method_index}(1, trial_index) = unc_stereo_individual{trial_index}{method_index}.x;
        unc_stereo_individual_all{method_index}(2, trial_index) = unc_stereo_individual{trial_index}{method_index}.y;
        unc_stereo_individual_all{method_index}(3, trial_index) = unc_stereo_individual{trial_index}{method_index}.z;
    end

    % uncertainties - combined
    for method_index = 1:num_combination_methods
        unc_stereo_combined_all{method_index}(1, trial_index) = unc_stereo_combined{trial_index}{method_index}.x;
        unc_stereo_combined_all{method_index}(2, trial_index) = unc_stereo_combined{trial_index}{method_index}.y;
        unc_stereo_combined_all{method_index}(3, trial_index) = unc_stereo_combined{trial_index}{method_index}.z;
    end
end

% ===============================
%% remove invalid measurements
% ===============================
fprintf('removing invalid measurements\n');

% errors
err_all_valid = cell(1, 4);
valid_indices = cell(1, 3);
err_all_valid{4} = [];
for i = 1:3
    % identify invalid measurements
    valid_indices{i} = find(abs(err_all(i, :)) > error_min & abs(err_all(i, :)) < error_max(i));
    err_all_valid{i} = err_all(i, valid_indices{i});
    err_all_valid{4} = [err_all_valid{4}, err_all_valid{i}];
end

% uncertainties - individual
unc_stereo_individual_all_valid = cell(num_uncertainty_methods, 4);
for method_index = 1:num_uncertainty_methods
    unc_stereo_individual_all_valid{method_index, 4} = [];
    for i = 1:3
        unc_stereo_individual_all_valid{method_index, i} = unc_stereo_individual_all{method_index}(i, valid_indices{i});
        unc_stereo_individual_all_valid{method_index, 4} = [unc_stereo_individual_all_valid{method_index, 4}, unc_stereo_individual_all_valid{method_index, i}];
    end    
end

% uncertainties - combined
unc_stereo_combined_all_valid = cell(num_combination_methods, 4);
for method_index = 1:num_combination_methods
    unc_stereo_combined_all_valid{method_index, 4} = [];
    for i = 1:3
        unc_stereo_combined_all_valid{method_index, i} = unc_stereo_combined_all{method_index}(i, valid_indices{i});
        unc_stereo_combined_all_valid{method_index, 4} = [unc_stereo_combined_all_valid{method_index, 4}, unc_stereo_combined_all_valid{method_index, i}];
    end
end

% ===============================
%% calculate pdf and cdf
% ===============================
fprintf('calculating pdf and cdf\n');

% errors
pdf_err = cell(1, 4);
cdf_err = cell(1, 4);
for i = 1:4
    pdf_err{i} = histcounts(abs(err_all_valid{i}), bins{i}, 'normalization', 'pdf');
    cdf_err{i} = histcounts(abs(err_all_valid{i}), bins{i}, 'normalization', 'cdf');
end

% uncertainties - individual
pdf_unc_individual = cell(num_uncertainty_methods, 4);
cdf_unc_individual = cell(num_uncertainty_methods, 4);
for method_index = 1:num_uncertainty_methods
    for i = 1:4
        pdf_unc_individual{method_index, i} = histcounts(unc_stereo_individual_all_valid{method_index, i}, bins{i}, 'normalization', 'pdf');
        cdf_unc_individual{method_index, i} = histcounts(unc_stereo_individual_all_valid{method_index, i}, bins{i}, 'normalization', 'cdf');
    end
end

% uncertainties - combined
pdf_unc_combined = cell(num_combination_methods, 4);
cdf_unc_combined = cell(num_combination_methods, 4);
for method_index = 1:num_combination_methods
    for i = 1:4
        pdf_unc_combined{method_index, i} = histcounts(unc_stereo_combined_all_valid{method_index, i}, bins{i}, 'normalization', 'pdf');
        cdf_unc_combined{method_index, i} = histcounts(unc_stereo_combined_all_valid{method_index, i}, bins{i}, 'normalization', 'cdf');
    end
end

% ===============================
%% calculate rms
% ===============================
fprintf('calculating rms\n');

% errors
rms_err = nans(1, 4);
for i = 1:4
    rms_err(i) = rms(err_all_valid{i});
end

% uncertainties - individual
rms_unc_individual = nans(num_uncertainty_methods, 4);
for method_index = 1:num_uncertainty_methods
    temp_valid_all = [];
    for i = 1:3
        % rms_unc_individual(method_index, i) = rms(unc_stereo_individual_all_valid{method_index, i});
        temp = unc_stereo_individual_all_valid{method_index, i};
        temp_valid = temp(temp < error_max(i));
        rms_unc_individual(method_index, i) = rms(temp_valid);
        temp_valid_all = [temp_valid_all, temp_valid];
    end
    rms_unc_individual(method_index, 4) = rms(temp_valid_all);
end

% uncertainties - combined
rms_unc_combined = nans(num_combination_methods, 4);
for method_index = 1:num_combination_methods
    temp_valid_all = [];
    for i = 1:3
        % rms_unc_combined(method_index, i) = rms(unc_stereo_combined_all_valid{method_index, i});
        temp = unc_stereo_combined_all_valid{method_index, i};
        temp_valid = temp(temp < error_max(i));
        rms_unc_combined(method_index, i) = rms(temp_valid);
        temp_valid_all = [temp_valid_all, temp_valid];
    end
    rms_unc_combined(method_index, 4) = rms(temp_valid_all);
end

% ==================
% save results
% ==================
fprintf('saving results to file\n');
filename = fullfile(current_read_directory, 'stereo_errors_uncertainties_statistics.mat');
save(filename, 'err_all', 'unc_stereo_individual_all', 'unc_stereo_combined_all', ...
                'valid_indices', 'err_all_valid', 'unc_stereo_individual_all_valid', 'unc_stereo_combined_all_valid', ...
                'pdf_err', 'cdf_err', 'pdf_unc_individual', 'cdf_unc_individual', 'pdf_unc_combined', 'cdf_unc_combined', ...
                'rms_err', 'rms_unc_individual', 'rms_unc_combined');





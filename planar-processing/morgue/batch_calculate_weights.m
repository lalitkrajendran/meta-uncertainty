clear
close all
clc

restoredefaultpath;
addpath(genpath('/scratch/shannon/c/aether/Projects/BOS/general-codes/matlab-codes/'));
addpath('../prana-uncertainty-average-dc-new-im-2-save-deform-cs/');
dbstop if error

%% read/write settings

% window resolution
window_resolution_array = [32, 64];
num_window_size = numel(window_resolution_array);

% pass number
pass_number = 4;

% directory containing images
top_image_directory = fullfile('/scratch/shannon/a/aether/Projects/PIV_PTV_Uncertainty_Quantification/Planar_uncertainty_work/Images/');
% directory containing files to be read
top_read_directory = fullfile('/scratch/shannon/a/aether/Projects/PIV_PTV_Uncertainty_Quantification/uncertainty-combination/analysis/results/', 'experiment-new');

% array of data set names
dataset_name_array = {'PivChal03B'; 'PivChal05B'; 'stagnation_flow'; 'Vortex_Ring'; 'Jetdata'};
num_datasets = numel(dataset_name_array);

% array of uncertainty methods
uncertainty_method_array = {'IM'; 'MC'; 'CS'};
num_uncertainty_methods = numel(uncertainty_method_array);

% array of uncertainty combination model to be used
% 'unweighted', 'global-weight-var', 'global-weight-std',
% 'local-weight-var', 'local-weight-std'
combination_methods = {'unwt'; 'var'; 'entropy'}; %'prob'};
num_combination_methods = numel(combination_methods);

% directory where results of this analysis are to be saved
top_write_directory = fullfile('/scratch/shannon/a/aether/Projects/PIV_PTV_Uncertainty_Quantification/uncertainty-combination/analysis/results/', 'monte-carlo', 'new-processing-all-datasets');
mkdir_c(top_write_directory);

%% analysis settings

% number of trials
num_trials = 1e3;
% minimum allowable error (pix.)
min_error_threshold = 1e-4;
% maximum allowable error (pix.)
max_error_threshold = 0.2;
% number of bins for histogram
num_bins = 30;
% bins for histograms
bins = linspace(min_error_threshold, max_error_threshold, num_bins);
% number of bins for the coarse histogram
num_bins_coarse = 8;
% bins for histogram of weights
bins_weights = linspace(0, 1, 25);

%% bootstrapping settings

% number of particles to remove
% num_particles_remove_array = [5, 10, 15];
percentage_particles_remove = 0.25;
% number of resampling trials
num_resampling_trials = 1e3;

% uncertainty flags for the various methods (0 = False, 1 = True)
uncertainty_flags.ppruncertainty = 0;
uncertainty_flags.miuncertainty = 0;
uncertainty_flags.mcuncertainty = 1;

%% particle identification settings

% intensity threshold for particle identification
intensity_threshold = 10;

% particle diameter (pix.)
d_p = 3;

% particle sizing settings
sizeprops = struct;
sizeprops.method = 'IWC';
sizeprops.p_area = 3; %0.5 * d_p^2;
sizeprops.sigma = 4;
sizeprops.errors = 0;
sizeprops.save_dir = [];

%% directory settings for this case

% directory to load results for this case
current_read_directory = fullfile(top_write_directory, ['trials=' num2str(num_trials, '%d')], ['ppr' num2str(percentage_particles_remove, '%.2f') '_nrs' num2str(num_resampling_trials)], 'new');

%% load results

load(fullfile(current_read_directory, 'monte_carlo_results.mat'));

%% loop through combination methods

unc_combined = cell(1, num_trials);
weights = cell(1, num_trials);

%% calculate weights and combined uncertainty

for trial_index = 1:num_trials
    fprintf('trial_index: %d\n', trial_index);
    %% extract results for this trial
    unc_current = unc_trials{trial_index};
    unc_resampling_current = unc_resampling{trial_index};

    % skip this trial if empty
    if isempty(unc_resampling_current)
        continue;
    end

    % initialize
    weights{trial_index} = cell(1, num_combination_methods);
    unc_combined{trial_index} = cell(1, num_combination_methods);
    %% unweighted
    combination_method_index = 1;
    
    % calculate weighting matrix
    weights_current = [1, 1, 1];
    weights{trial_index}{combination_method_index}.x = real(weights_current/sum(weights_current));

    weights_current = [1, 1, 1];
    weights{trial_index}{combination_method_index}.y = real(weights_current/sum(weights_current));
    
    %% variance covariance
    combination_method_index = 2;
    
    % calculate covariance matrix
    covariance_matrix_x = cov([unc_resampling_current.imx', unc_resampling_current.mcx', unc_resampling_current.csx']);
    covariance_matrix_y = cov([unc_resampling_current.imy', unc_resampling_current.mcy', unc_resampling_current.csy']);

    % calculate weighting matrix
    weights_current = [sum(1./covariance_matrix_x(:, 1)), sum(1./covariance_matrix_x(:, 2)), sum(1./covariance_matrix_x(:, 3))];
%     weights_current = [sum(1./covariance_matrix_x(1, 1)), sum(1./covariance_matrix_x(2, 2)), sum(1./covariance_matrix_x(3, 3))];
    weights{trial_index}{combination_method_index}.x = real(weights_current/sum(weights_current(:)));

    weights_current = [sum(1./covariance_matrix_y(:, 1)), sum(1./covariance_matrix_y(:, 2)), sum(1./covariance_matrix_y(:, 3))];
%     weights_current = [sum(1./covariance_matrix_y(1, 1)), sum(1./covariance_matrix_y(2, 2)), sum(1./covariance_matrix_y(3, 3))];
    weights{trial_index}{combination_method_index}.y = real(weights_current/sum(weights_current(:)));

    %% entropy    
    combination_method_index = 3;

    entropy = nans(3, 2);
    % -------------------------
    % calculate entropies
    % -------------------------
    entropy(1, 1) = calculate_shannon_entropy(unc_resampling_current.imx, bins);
    entropy(1, 2) = calculate_shannon_entropy(unc_resampling_current.imy, bins);
    entropy(2, 1) = calculate_shannon_entropy(unc_resampling_current.mcx, bins);
    entropy(2, 2) = calculate_shannon_entropy(unc_resampling_current.mcy, bins);
    entropy(3, 1) = calculate_shannon_entropy(unc_resampling_current.csx, bins);
    entropy(3, 2) = calculate_shannon_entropy(unc_resampling_current.csy, bins);
    
    % -------------------------
    % calculate weights
    % -------------------------
    weights_current = 1./entropy(:, 1);
    weights{trial_index}{combination_method_index}.x = real(weights_current/sum(weights_current));

    weights_current = 1./entropy(:, 2);
    weights{trial_index}{combination_method_index}.y = real(weights_current/sum(weights_current));

    %% probability
%     combination_method_index = 4;
% 
%     probability = nans(3, 2);
%     % -------------------------
%     % calculate entropies
%     % -------------------------
%     probability(1, 1) = calculate_probability(unc_resampling_current.imx, unc_current.imx, bins); 
%     probability(1, 2) = calculate_probability(unc_resampling_current.imy, unc_current.imy, bins); 
%     probability(2, 1) = calculate_probability(unc_resampling_current.mcx, unc_current.mcx, bins); 
%     probability(2, 2) = calculate_probability(unc_resampling_current.mcy, unc_current.mcy, bins); 
%     probability(3, 1) = calculate_probability(unc_resampling_current.csx, unc_current.csx, bins); 
%     probability(3, 2) = calculate_probability(unc_resampling_current.csy, unc_current.csy, bins); 
%     
%     % -------------------------
%     % calculate weights
%     % -------------------------
%     weights_current = probability(:, 1); 
%     weights{trial_index}{combination_method_index}.x = real(weights_current/sum(weights_current));
% 
%     weights_current = probability(:, 2);
%     weights{trial_index}{combination_method_index}.y = real(weights_current/sum(weights_current));
    
    %% calculate the combined uncertainty
    
    for combination_method_index = 1:num_combination_methods
        unc_combined{trial_index}{combination_method_index} = struct;
        unc_combined{trial_index}{combination_method_index}.x = dot([unc_current.imx, unc_current.mcx, unc_current.csx], weights{trial_index}{combination_method_index}.x);
        unc_combined{trial_index}{combination_method_index}.y = dot([unc_current.imy, unc_current.mcy, unc_current.csy], weights{trial_index}{combination_method_index}.y);
    end
end

%% save calculations to file

fprintf('saving results to file\n');
filename = fullfile(current_read_directory, 'combined_uncertainties.mat');
save(filename, 'weights', 'unc_combined');

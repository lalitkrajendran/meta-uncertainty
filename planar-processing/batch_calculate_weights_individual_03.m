% This script performs resampling calculations for a range of particle removal
% percentage and plots the change in uncertainties for different methods

clear
close all
clc

restoredefaultpath;
addpath(genpath('/scratch/shannon/c/aether/Projects/BOS/general-codes/matlab-codes/'));
addpath('../prana/')
addpath('../general-codes/')
addpath('../histogram_distance/')
setup_default_settings;

% ============================
%% read/write settings
% ============================
% window resolution
window_resolution_array = [64, 32]; 
num_window_resolution = numel(window_resolution_array);
% pass number
pass_number = 4;
% directory containing images
top_image_directory = fullfile('/scratch/shannon/a/aether/Projects/PIV_PTV_Uncertainty_Quantification/Planar_uncertainty_work/Images/');
% directory containing files to be read
top_read_directory = fullfile('/scratch/shannon/a/aether/Projects/PIV_PTV_Uncertainty_Quantification/uncertainty-combination/analysis/results/planar-dataset', 'experiment-new');
% array of data set names
dataset_name_array = {'PivChal03B'; 'PivChal05B'; 'stagnation_flow'; 'Vortex_Ring'; 'Jetdata'};
num_datasets = numel(dataset_name_array);
% array of uncertainty methods
individual_method_array = {'IM'; 'MC'; 'CS'};
num_individual_methods = numel(individual_method_array);
% array of snr metrics
snr_metric_array = {'PPR'; 'MI'};
num_snr_methods = numel(snr_metric_array);
% components
components = {'x'; 'y'};
num_components = numel(components);
% directory where results of this analysis are to be saved
top_write_directory = fullfile('/scratch/shannon/a/aether/Projects/PIV_PTV_Uncertainty_Quantification/uncertainty-combination/analysis/results/planar-dataset', 'monte-carlo', 'individual-datasets');
mkdir_c(top_write_directory);
% number of trials
num_trials = 1e3;

% ============================
%% resampling settings
% ============================
% resampling method names
resampling_method_names = {'remove-paired'; 'remove-random'; 'add-random'};
% resampling_method_names = {'add-random'};
num_resampling_methods = numel(resampling_method_names);
% number of particles to remove
percentage_particles_remove_array = 0.05:0.05:0.25;
num_ppr = numel(percentage_particles_remove_array);
% number of resampling trials
num_resampling_trials = 1e2; %3;
% display_resampled_images? (true/false)
display_resampled_images = 0;
% method names for plotting
resampling_method_names_plot = {'Removing Paired Particles'; 'Removing Random Particles'; 'Adding Random Particles'};
resampling_method_names_plot_short = {'R-P', 'R-R', 'A-R'};
% % case name
% resampling_case_name = ['ppr' num2str(percentage_particles_remove, '%.2f') '_nrs' num2str(num_resampling_trials) '-prewindow-randpart-snr-newproc'];

% ========================
%% weight settings
% ========================
% number of bins for weights
num_bins_w = 30;
% bins for weights
bins_w = linspace(0, 1, num_bins_w);
% metric type ('rms', 'iqr', 'std')
metric_name = 'iqr';

% ========================
%% statistics settings
% ========================
min_error_threshold = 1e-3;
max_error_threshold = 0.2;
num_bins = 30;
bins = linspace(-max_error_threshold, max_error_threshold, num_bins*2);
bins_abs = linspace(1e-3, max_error_threshold, num_bins);
% methods to calculate histogram distances
histogram_distance_method = 'total_variation_distance';
% categories
histogram_distance_categories = categorical({individual_method_array{:}, resampling_method_names_plot_short{:}});
% sort in the desired order
histogram_distance_categories = reordercats(histogram_distance_categories, ...
                                {individual_method_array{:}, resampling_method_names_plot_short{:}});
% error model ('gaussian' or 'lognormal')
% error_models = {'gaussian'; 'lognormal'; 'exp'};
% error_models = {'rayleigh'; 'lognormal'; 'exp'};
error_models = {'gaussian'};
num_models = numel(error_models);

% ============================
%% pre-load all dataset
% ============================
fprintf('Loading all datasets into memory\n');

results_all = cell(num_window_resolution, num_datasets);
jobfile_all = cell(num_window_resolution, num_datasets);
files_im1 = cell(num_window_resolution, num_datasets);
files_im2 = cell(num_window_resolution, num_datasets);

% loop through window resolutions
for window_resolution_index = 1:num_window_resolution
    % loop through datasets
    for dataset_index = 1:num_datasets
        % name of the data set
        dataset_name = dataset_name_array{dataset_index};        
        fprintf('Dataset: %s\n', dataset_name);

        % ============================
        %% load data
        % ============================
        % results directory for current data set
        current_results_directory = fullfile(top_read_directory, dataset_name, ['WS' num2str(window_resolution_index)]);

        % directory containing vectors
        vectors_directory = fullfile(current_results_directory, 'vectors-new');

        % load results for vectors, errors and uncertainties
        results_all{window_resolution_index, dataset_index} = load_directory_data(vectors_directory, ['*pass' num2str(pass_number) '*.mat']);

        % load jobfile
        jobfile_all{window_resolution_index, dataset_index} = load(fullfile(current_results_directory, 'jobfile.mat'));
        
        % ============================
        %% load listing of deformed images
        % ============================        
        % directory containing deformed images for current data set
        deformed_images_directory = fullfile(vectors_directory, 'imDeform');
        
        % get list of im1 files in the directory
        files_im1{window_resolution_index, dataset_index} = get_directory_listing(deformed_images_directory, 'PIV*im1d*.mat');
        % get list of im2 files in the directory
        files_im2{window_resolution_index, dataset_index} = get_directory_listing(deformed_images_directory, 'PIV*im2d*.mat');        
    end
end

% ============================
%% load errors for all datasets
% ============================
fprintf('Loading all errors into memory\n');
errors_all = cell(num_window_resolution, num_datasets);

% loop through window resolutions
for window_resolution_index = 1:num_window_resolution
    % loop through datasets
    for dataset_index = 1:num_datasets
        % name of the data set        
        dataset_name = dataset_name_array{dataset_index};
        
        fprintf('Dataset: %s\n', dataset_name);
        
        %% Load data

        % results directory for current data set
        current_results_directory = fullfile(top_read_directory, dataset_name, ['WS' num2str(window_resolution_index)]);

        % load results for vectors, errors and uncertainties
        errors_all{window_resolution_index, dataset_index} = load(fullfile(current_results_directory, 'errors-new.mat'));        
    end
end

fprintf('running monte-carlo\n');
% set seed for random number generator
rng(0);

X_all = cell(num_window_resolution, num_datasets);
Y_all = cell(num_window_resolution, num_datasets);
unc_mean_all = cell(num_window_resolution, num_datasets);
d_unc_mean_all = cell(num_window_resolution, num_datasets);
snr_mean_all = cell(num_window_resolution, num_datasets);

% ============================
%% loop through window resolutions
% ============================
for window_resolution_index = 1:num_window_resolution
    fprintf('WS: %d\n', window_resolution_index);

    % ============================
    %% loop through datasets
    % ============================
    for dataset_index = 1:num_datasets
        % dataset name
        dataset_name = dataset_name_array{dataset_index};
        fprintf('dataset: %s\n', dataset_name);

        % ================================================
        % load results
        % ================================================
        % directory to save results for this case
        current_read_directory = fullfile(top_write_directory, dataset_name, ['WS' num2str(window_resolution_index)], 'weights-rms-change-study');
        % file name
        filename = ['resampling-statistics-nt=' num2str(num_trials, '%d') 'nr='  num2str(num_resampling_trials, '%d') '.mat'];
        % load results to file
        load(fullfile(current_read_directory, filename));
        
        % --------------------------
        % initialize 
        % --------------------------            
        unc_ratio = cell(num_resampling_methods, num_trials);
        unc_ratio_fit = cell(1, num_resampling_methods, num_trials);
        unc_ratio_rate = cell(1, num_resampling_methods, num_trials);
        wt_temp = cell(num_resampling_methods, num_trials);        

        % ==========================
        % loop through trials
        % ==========================
        for trial_index = 1:num_trials
            if rem(trial_index, 100) == 0
                fprintf('trial_index: %d\n', trial_index);
            end
            % ==========================
            % loop through resampling methods
            % ==========================
            for resampling_method_index = 1:num_resampling_methods          
                % calculate weight
                [wt_temp{resampling_method_index, trial_index}, unc_ratio{resampling_method_index, trial_index}, unc_ratio_fit{resampling_method_index, trial_index}] = calculate_weights_from_resampled_uncertainty_ratios(unc_sub_trials{trial_index}, {unc_resampling_trials{trial_index, :, resampling_method_index}}, ...
                                                                metric_name, individual_method_array, percentage_particles_remove_array, components, num_resampling_trials, 'struct');
            end    
        end

        % ==========================
        % re-arrange weights to be consistent with rest of the code
        % ==========================        
        wt_trials = cell(1, num_resampling_methods);
        for resampling_method_index = 1:num_resampling_methods
            wt_trials{resampling_method_index} = struct;
            for component_index = 1:num_components
                component_name = components{component_index};
                wt_trials{resampling_method_index}.(component_name) = nans(num_trials, num_individual_methods);
                for trial_index = 1:num_trials
                    wt_trials{resampling_method_index}.(component_name)(trial_index, :) = wt_temp{resampling_method_index, trial_index}.(component_name);
                end
            end          
        end

        % ==========================
        % extract errors
        % ==========================
        errors = errors_all{window_resolution_index, dataset_index};
        % calculate linear indices for the grid points across trials
        idx = sub2ind(size(errors.err_U), r_trials, c_trials, snapshot_trials);

        % extract errors
        err.x = errors.err_U(idx)';
        err.y = errors.err_V(idx)';

        % ==========================
        % extract individual uncertainty
        % ==========================
        unc_indiv.x = nans(num_trials, num_individual_methods);
        unc_indiv.y = nans(num_trials, num_individual_methods);
        for trial_index = 1:num_trials
            results = results_all{window_resolution_index, dataset_index}{snapshot_trials(trial_index)};
            unc_current = extract_planar_uncertainties(results.uncertainty2D, r_trials(trial_index), c_trials(trial_index));
            for individual_method_index = 1:num_individual_methods
                method_name = lower(individual_method_array{individual_method_index});
                unc_indiv.x(trial_index, individual_method_index) = unc_current.([method_name 'x']);
                unc_indiv.y(trial_index, individual_method_index) = unc_current.([method_name 'y']);
            end
        end

        % ==========================
        % calculate combined uncertainty
        % ==========================
        unc_comb.x = nans(num_trials, num_resampling_methods);
        unc_comb.y = nans(num_trials, num_resampling_methods);
        for resampling_method_index = 1:num_resampling_methods
            unc_comb.x(:, resampling_method_index) = sum(unc_indiv.x .* wt_trials{resampling_method_index}.x, 2); 
            unc_comb.y(:, resampling_method_index) = sum(unc_indiv.y .* wt_trials{resampling_method_index}.y, 2); 
        end                
        
        % ==========================
        % save results
        % ==========================
        % directory to save results for this case
        current_read_directory = fullfile(top_write_directory, dataset_name, ['WS' num2str(window_resolution_index)], 'weights-rms-change-study');
        % file name
        filename = ['resampling-statistics-nt=' num2str(num_trials, '%d') 'nr='  num2str(num_resampling_trials, '%d') '-weights-unc-' metric_name '-03.mat'];
        % load results to file
        save(fullfile(current_read_directory, filename), 'unc_ratio', 'unc_ratio_fit', 'err', 'unc_indiv', 'unc_comb', 'wt_trials'); 
    end
end

shut_down_parpool();

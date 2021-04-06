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
        current_read_directory = fullfile(top_write_directory, dataset_name, ['WS' num2str(window_resolution_index)], 'weights-rms-change-study-profile');
        % file name
        filename = ['resampling-statistics-nr='  num2str(num_resampling_trials, '%d') '.mat'];
        % load results to file
        load(fullfile(current_read_directory, filename));

        % load errors
        errors = errors_all{window_resolution_index, dataset_index};
        
        % calculate number of snapshots available for the current dataset and
        % window size
        num_snapshots = numel(results_all{window_resolution_index, dataset_index});
        if num_snapshots > size(errors.err_U, 3)
            num_snapshots = size(errors.err_U, 3);
        end

        % calculate number of grid points
        num_grid_points_prof = numel(r_prof);

        % --------------------------
        % initialize 
        % --------------------------            
        unc_ratio = cell(num_resampling_methods, num_grid_points_prof);
        unc_ratio_fit = cell(1, num_resampling_methods, num_grid_points_prof);
        unc_ratio_rate = cell(1, num_resampling_methods, num_grid_points_prof);
        wt_temp = cell(num_resampling_methods, num_grid_points_prof);        
        
        % ==========================
        % loop through resampling methods
        % ==========================
        for resampling_method_index = 1:num_resampling_methods
            fprintf('resampling method: %d of %d\n', resampling_method_index, num_resampling_methods);
            % ==========================
            % loop through snapshots
            % ==========================
            for snapshot_index = 1:num_snapshots
                fprintf('snapshot: %d of %d\n', snapshot_index, num_snapshots);
                % ================================================
                %% loop through grid points in the profile
                % ================================================
                parfor grid_point_index = 1:num_grid_points_prof                
                    % fprintf('grid point: %d of %d\n', grid_point_index, num_grid_points_prof);                
                    % calculate weight
                    [wt, unc_rat, unc_rat_fit] = calculate_weights_from_resampled_uncertainty_ratios(unc_sub_trials{snapshot_index, grid_point_index}, ...
                                                                {unc_resampling_trials{snapshot_index, grid_point_index, :, resampling_method_index}}, ...
                                                                    metric_name, individual_method_array, percentage_particles_remove_array, ...
                                                                    components, num_resampling_trials, 'struct');

                    wt_temp{resampling_method_index, snapshot_index, grid_point_index} = wt;
                    unc_ratio{resampling_method_index, snapshot_index, grid_point_index} = unc_rat;
                    unc_ratio_fit{resampling_method_index, snapshot_index, grid_point_index} = unc_rat_fit;
                    % [wt_temp{resampling_method_index, grid_point_index}, unc_ratio{resampling_method_index, grid_point_index}, unc_ratio_fit{resampling_method_index, grid_point_index}] = calculate_weights_from_resampled_uncertainty_ratios(unc_sub_trials{trial_index}, {unc_resampling_trials{trial_index, :, resampling_method_index}}, ...
                    %                                                 metric_name, individual_method_array, percentage_particles_remove_array, components, num_resampling_trials, 'struct');
                end     
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
                wt_trials{resampling_method_index}.(component_name) = nans(num_snapshots, num_grid_points_prof, num_individual_methods);
                for snapshot_index = 1:num_snapshots
                    for grid_point_index = 1:num_grid_points_prof
                        wt_trials{resampling_method_index}.(component_name)(snapshot_index, grid_point_index, :) = wt_temp{resampling_method_index, snapshot_index, grid_point_index}.(component_name);
                    end
                end
            end          
        end        

        % ==========================
        % extract errors
        % ==========================
        err = struct;
        err.x = nans(num_snapshots, num_grid_points_prof);
        err.y = nans(num_snapshots, num_grid_points_prof);
        
        for snapshot_index = 1:num_snapshots
            % calculate linear indices for the grid points across trials
            idx = sub2ind(size(errors.err_U), r_prof, c_prof, snapshot_index*ones(size(r_prof)));
            
            err.x(snapshot_index, :) = errors.err_U(idx);
            err.y(snapshot_index, :) = errors.err_V(idx);
        end

        % ==========================
        % extract individual uncertainty
        % ==========================
        unc_indiv = struct;
        for component_index = 1:num_components
            component_name = components{component_index};
            unc_indiv.(component_name) = nans(num_snapshots, num_grid_points_prof, num_individual_methods);
            for snapshot_index = 1:num_snapshots
                for grid_point_index = 1:num_grid_points_prof
                    unc_current = unc_orig{snapshot_index, grid_point_index};
                    for individual_method_index = 1:num_individual_methods
                        method_name = lower(individual_method_array{individual_method_index});
                        unc_indiv.(component_name)(snapshot_index, grid_point_index, individual_method_index) = unc_current.([method_name component_name]);
                    end    
                end
            end
        end

        % ==========================
        % calculate combined uncertainty
        % ==========================
        unc_comb = struct;
        for component_index = 1:num_components
            component_name = components{component_index};
            unc_comb.(component_name) = nans(num_snapshots, num_grid_points_prof, num_resampling_methods);
            for resampling_method_index = 1:num_resampling_methods
                unc_comb.(component_name)(:, :, resampling_method_index) = sum(unc_indiv.(component_name) .* wt_trials{resampling_method_index}.(component_name), 3); 
            end
        end
        % ==========================
        % save results
        % ==========================
        % directory to save results for this case
        current_read_directory = fullfile(top_write_directory, dataset_name, ['WS' num2str(window_resolution_index)], 'weights-rms-change-study-profile');
        % file name
        filename = ['resampling-statistics-nr='  num2str(num_resampling_trials, '%d') '-weights-unc-' metric_name '-03.mat'];
        % load results to file
        save(fullfile(current_read_directory, filename), 'unc_ratio', 'unc_ratio_fit', 'err', 'unc_indiv', 'unc_comb', 'wt_trials'); 
    end
end

shut_down_parpool();

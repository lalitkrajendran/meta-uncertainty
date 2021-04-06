% This script calculates weights for uncertainty distributions using
% different combination methods.

clear
close all
clc

restoredefaultpath;
addpath(genpath('/scratch/shannon/c/aether/Projects/BOS/general-codes/matlab-codes/'));
% addpath('../prana-uncertainty-average-dc-new-im-2-save-deform-cs/');
addpath('../prana/');
addpath('../general-codes') 
addpath('../CompPD/')
setup_default_settings;

% ===============================
%% read/write settings
% ===============================
% window resolution
window_resolution_array = [32, 64];
num_window_resolution = numel(window_resolution_array);

% pass number
pass_number = 4;

% directory containing images
top_image_directory = fullfile('/scratch/shannon/a/aether/Projects/PIV_PTV_Uncertainty_Quantification/Planar_uncertainty_work/Images/');
% directory containing files to be read
top_read_directory = fullfile('/scratch/shannon/a/aether/Projects/PIV_PTV_Uncertainty_Quantification/uncertainty-combination/analysis/results/planar-dataset', 'experiment-new');

% array of data set names
dataset_name_array = {'PivChal03B'}; %'PivChal05B'; 'stagnation_flow'; 'Vortex_Ring'; 'Jetdata'};
num_datasets = numel(dataset_name_array);

% array of uncertainty methods
individual_method_array = {'IM'; 'MC'; 'CS'};
num_individual_methods = numel(individual_method_array);

% array of snr metrics
snr_metric_array = {'PPR'; 'MI'};
num_snr_methods = numel(snr_metric_array);

% array of uncertainty combination model to be used
% 'unweighted', 'global-weight-var', 'global-weight-std',
% 'local-weight-var', 'local-weight-std'
% combination_methods = {'unwt'; 'pd-var'; 'entropy'}; %'prob'};
combination_methods = {'unwt'; 'var'; 'entropy'}; %'prob'};
num_combination_methods = numel(combination_methods);

% component names
component_names = {'x'; 'y'};
num_components = numel(component_names);

% directory where results of this analysis are to be saved
top_write_directory = fullfile('/scratch/shannon/a/aether/Projects/PIV_PTV_Uncertainty_Quantification/uncertainty-combination/analysis/results/planar-dataset', ...
                                'monte-carlo', 'individual-datasets');
mkdir_c(top_write_directory);

% ===============================
%% analysis settings
% ===============================
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

% ===============================
%% resampling settings
% ===============================
% number of particles to remove
% num_particles_remove_array = [5, 10, 15];
percentage_particles_remove = 0.25;
% number of resampling trials
num_resampling_trials = 1e3;
% case name
resampling_case_name = ['ppr' num2str(percentage_particles_remove, '%.2f') '_nrs' num2str(num_resampling_trials) '-prewindow-randpart-snr-newproc'];

% ===============================
%% loop through window resolutions
% ===============================
for window_resolution_index = 1:num_window_resolution
    fprintf('WS: %d\n', window_resolution_index);
    
    % ===============================
    %% loop through datasets
    % ===============================
    for dataset_index = 1:num_datasets
        fprintf('dataset: %s\n', dataset_name_array{dataset_index});
        % directory to load results for this case
        current_read_directory = fullfile(top_write_directory, dataset_name_array{dataset_index}, ['WS' num2str(window_resolution_index)], ...
        ['trials=' num2str(num_trials, '%d')], resampling_case_name);

        % load results
        load(fullfile(current_read_directory, 'monte_carlo_results.mat'));

        % ================================================
        % aggregate resampling results
        % ================================================
        % uncertainties        
        fprintf('aggregating resampled uncertainties\n');
        unc_resampled = cell(1, num_trials);
        for trial_index = 1:num_trials
            unc_resampled{trial_index} = cell(1, num_individual_methods);
            for individual_method_index = 1:num_individual_methods
                unc_resampled{trial_index}{individual_method_index} = struct;
                unc_resampled{trial_index}{individual_method_index}.name = individual_method_array{individual_method_index};
                for component_index = 1:num_components
                    unc_resampled{trial_index}{individual_method_index}.(component_names{component_index}) = ...
                                        unc_resampling{trial_index}.([lower(individual_method_array{individual_method_index}) component_names{component_index}]);
                end
            end
        end

        % snr metric
        fprintf('aggregating resampled snr metric\n');
        snr_resampled = cell(1, num_trials);
        for trial_index = 1:num_trials
            snr_resampled{trial_index} = cell(1, num_snr_methods);
            for snr_method_index = 1:num_snr_methods
                snr_resampled{trial_index}{snr_method_index} = repmat(snr_metric_resampling{trial_index}.(snr_metric_array{snr_method_index}), 1, num_components);
            end
        end
        
        % ================================================
        % remove invalid resampled uncertainties
        % ================================================
        fprintf('removing invalid resampled uncertainties\n');
        for trial_index = 1:num_trials
            for individual_method_index = 1:num_individual_methods
                for component_index = 1:num_components                    
                    % extract current uncertainties
                    unc_current = unc_resampled{trial_index}{individual_method_index}.(component_names{component_index});                    
                    % nan invalid measurements
                    unc_current = nan_invalid_measurements(unc_current, min_error_threshold, max_error_threshold);
                    % unc_current = remove_invalid_measurements(unc_current, min_error_threshold, max_error_threshold);
                    % update uncertainties
                    unc_resampled{trial_index}{individual_method_index}.(component_names{component_index}) = unc_current;
                end
            end
        end

        % ===============================
        %% calculate weights and combined uncertainty
        % ===============================
        fprintf('calculating weights and combined uncertainty\n');
        unc_combined = cell(1, num_trials);
        weights = cell(1, num_trials);
        parfor trial_index = 1:num_trials
            fprintf('trial_index: %d\n', trial_index);
            
            % ============================
            % calculate weights
            % ============================
            weights{trial_index} = calculate_weights_from_resampling_general(unc_resampled{trial_index}, component_names, {bins; bins}, combination_methods);

            % ============================
            %% calculate the combined uncertainty
            % ============================            
            for combination_method_index = 1:num_combination_methods
                unc_combined{trial_index}{combination_method_index} = struct;
                unc_combined{trial_index}{combination_method_index}.name = combination_methods{combination_method_index};
                for component_index = 1:num_components
                    component_current = component_names{component_index};
                    unc_combined{trial_index}{combination_method_index}.(component_current) = ...
                                        [unc_trials{trial_index}.(['im' component_current]), ...
                                                unc_trials{trial_index}.(['mc' component_current]), ...
                                                unc_trials{trial_index}.(['cs' component_current])] * ...
                                                 weights{trial_index}{combination_method_index}.(component_current);
                end
            end
        end

        % ============================
        %% save calculations to file
        % ============================
        fprintf('saving results to file\n');
        filename = fullfile(current_read_directory, 'combined_uncertainties.mat');
        save(filename, 'weights', 'unc_combined', 'unc_resampled', 'snr_resampled');

        return;
    end
end

shut_down_parpool();
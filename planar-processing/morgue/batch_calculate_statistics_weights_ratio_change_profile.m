% This script performs resampling calculations for a range of particle removal
% percentage and plots the change in uncertainties for different methods

clear
close all
clc

restoredefaultpath;
addpath(genpath('/scratch/shannon/c/aether/Projects/BOS/general-codes/matlab-codes/'));
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
% dataset_name_array = {'Vortex_Ring'; 'Jetdata'};
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

% ============================
%% particle identification settings
% ============================
% intensity threshold for particle identification
intensity_threshold = 100;
% particle diameter (pix.)
d_p = 3;
% particle sizing settings
sizeprops = struct;
sizeprops.method = 'IWC';
sizeprops.p_area = 2; % this will be modified for each image
sizeprops.sigma = 4;
sizeprops.errors = 0;
sizeprops.save_dir = [];
% display id results?
display_id_results = 0;

% ========================
%% weight settings
% ========================
% number of bins for weights
num_bins_w = 30;
% bins for weights
bins_w = linspace(0, 1, num_bins_w);

% ========================
%% statistics settings
% ========================
min_error_threshold = 0;
max_error_threshold = 1; %00; %0.2;
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

% metric type ('rms', 'iqr', 'std')
metric_name = 'iqr';

% ==========================
% initialize variables
% ==========================
err_all_consolidated = [];
unc_indiv_all_consolidated = [];
unc_comb_all_consolidated = [];
valid_points_consolidated = [];
num_valid_points_datasets = [];
err_est_indiv_consolidated = [];
err_est_comb_consolidated = [];
wt_consolidated = cell(1, num_resampling_methods);

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
        % file name for weights and combined uncertainties
        filename = ['resampling-statistics-nr='  num2str(num_resampling_trials, '%d') '-weights-unc-' metric_name '-03.mat'];
        % load results to file
        load(fullfile(current_read_directory, filename));

        num_snapshots = size(err.x, 1);
        num_grid_points_prof = size(err.x, 2);
        num_points = numel(err.x);
    
        % ==========================
        % calculate rms across snapshots
        % ==========================
        err_rms = struct;
        unc_indiv_rms = struct;
        unc_comb_rms = struct;
        for component_index = 1:num_components
            % component name
            component_name = components{component_index};
            
            % error
            temp = err.(component_name);
            temp(abs(temp) > max_error_threshold) = NaN;
            err_rms.(component_name) = nanrms(temp, 1);

            % uncertainty, individual            
            temp = unc_indiv.(component_name);
            temp(abs(temp) > max_error_threshold) = NaN;
            unc_indiv_rms.(component_name) = squeeze(nanrms(temp, 1));

            % uncertainty, combined
            temp = unc_comb.(component_name);
            temp(abs(temp) > max_error_threshold) = NaN;
            unc_comb_rms.(component_name) = squeeze(nanrms(temp, 1));            
        end

        % ==========================
        % calculate median weight across snapshots
        % ==========================
        wt_med = cell(1, num_resampling_methods);
        for resampling_method_index = 1:num_resampling_methods
            for individual_method_index = 1:num_individual_methods
                wt_med{resampling_method_index} = struct;
                for component_index = 1:num_components
                    % component name
                    component_name = components{component_index};
                    % calculate median weight
                    wt_med{resampling_method_index}.(component_name) = squeeze(median(wt_trials{resampling_method_index}.(component_name)));
                end
            end
        end

        % ==========================
        % save results
        % ==========================
        % file name
        filename = ['resampling-statistics-nr='  num2str(num_resampling_trials, '%d') '-statistics-' metric_name '-new-thresh=' num2str(max_error_threshold, '%.2f') 'pix.mat'];
        % load results to file
        save(fullfile(current_read_directory, filename), 'wt_med', 'err_rms', 'unc_indiv_rms', 'unc_comb_rms');

    end
end
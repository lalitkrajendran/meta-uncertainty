% This script performs resampling calculations for a range of particle removal
% percentage and plots the change in uncertainties for different methods

clear
close all
clc

restoredefaultpath;
addpath(genpath('/scratch/shannon/c/aether/Projects/BOS/general-codes/matlab-codes/'));
addpath('../general-codes/')

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
% dataset_name_array_plot = {'03B'; '05B'; 'SF'; 'VR'; 'Jet'};
dataset_name_array_plot = {'TBL'; 'LSB'; 'SF'; 'VR'; 'Jet'};
% array of uncertainty methods
individual_method_array = {'IM'; 'MC'; 'CS'};
num_individual_methods = numel(individual_method_array);
% array of snr metrics
snr_metric_array = {'PPR'; 'MI'};
num_snr_methods = numel(snr_metric_array);
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
% resampling_method_names_plot_short = {'R-P', 'R-R', 'A-R'};
% resampling_method_names_plot_short = {'A-R'};
resampling_method_names_plot_short = {'Comb'};
resampling_method_index_plot = 3;

% % case name
% resampling_case_name = ['ppr' num2str(percentage_particles_remove, '%.2f') '_nrs' num2str(num_resampling_trials) '-prewindow-randpart-snr-newproc'];

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
max_error_threshold = 1.0;
max_error_plot = 0.2;
num_bins = 30;
bins = linspace(-max_error_plot, max_error_plot, num_bins*2);
bins_abs = linspace(1e-3, max_error_plot, num_bins);
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

% ========================
%% plot settings
% ========================
% display figures?
display_figures = 1;
% user screen resolution
user_screen_resolution = 113;
% save figure? (true/false)
save_figures = 1;
% symbols
symbols = {'o', '^', 'v'};
% colors
colors = lines(3);
% line symbols
line_symbols = {'-'; ':'; '-.'};
% range of displacements to be displayed in the contour plots
displacement_color_min = [0, 0, -0.25, -2, 0];
displacement_color_max = [15, 5, 0.25, 2, 5]; 

results = cell(num_window_resolution, num_datasets);
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

        % --------------------------
        % load statistics
        % --------------------------
        % file name
        % filename = ['resampling-statistics-nt=' num2str(num_trials, '%d') 'nr='  num2str(num_resampling_trials, '%d') '-statistics-' metric_name '-new.mat'];
        filename = ['resampling-statistics-nt=' num2str(num_trials, '%d') 'nr='  num2str(num_resampling_trials, '%d') '-statistics-' metric_name '-new-thresh=' num2str(max_error_threshold, '%.2f') 'pix.mat'];
        % load results to file
        results{window_resolution_index, dataset_index} = load(fullfile(current_read_directory, filename));
        % results{window_resolution_index+1, dataset_index} = results{window_resolution_index, dataset_index};
    end
end

% directory to save results for this case
write_directory = fullfile(top_write_directory, 'weights-rms-change-study-consolidated');
% create directory to save figures
if save_figures
    % figure_save_directory = fullfile(write_directory, ['figures-' metric_name]);
    figure_save_directory = fullfile(write_directory, ['figures-' metric_name '-thresh=' num2str(max_error_threshold, '%.2f') 'pix']);
    mkdir_c(figure_save_directory);
end

% ==========================
% plot violins
% ==========================
figure
% plot_consolidated_violins_v4(results, dataset_name_array_plot, individual_method_array, resampling_method_names_plot_short, max_error_threshold, user_screen_resolution);                   
plot_consolidated_violins_v5(results, dataset_name_array_plot, individual_method_array, resampling_method_names_plot_short, resampling_method_index_plot, max_error_plot, user_screen_resolution);                   
drawnow();

if save_figures
    save_figure_to_png_svg_fig(figure_save_directory, 'violin-panel', [1, 0, 0]);
end

return;
% ===================================
%% quantile-quantile plot 
% ===================================
error_model = 'gaussian'; 
figure
% plot_consolidated_qq_v2(results, dataset_name_array_plot, individual_method_array, resampling_method_names_plot_short, max_error_threshold, symbols, user_screen_resolution)
plot_consolidated_qq_v3(results, dataset_name_array_plot, individual_method_array, resampling_method_names_plot_short, resampling_method_index_plot, max_error_plot, symbols, user_screen_resolution)
drawnow();

% save figures
if save_figures
    save_figure_to_png_svg_fig(figure_save_directory, ['qq-error-true-vs-est-' error_model '-panel-v2'], [1, 0, 0]);
end

% % ===================================
% % total variation distance 
% % ===================================
% figure
% % plot_consolidated_distance_02(results, dataset_name_array_plot, individual_method_array, resampling_method_names_plot_short)
% plot_consolidated_distance_03(results, dataset_name_array_plot, individual_method_array, resampling_method_names_plot_short, resampling_method_index_plot)
% drawnow();
% % save figure
% if save_figures
%     save_figure_to_png_svg_fig(figure_save_directory, ['histogram-distance-' error_model '-panel'], [1, 0, 0]);
% else
%     pause(0.1);
% end    


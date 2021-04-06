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
dataset_name_array_plot = {'TBL'; 'LSB'; 'SF'; 'VR'; 'Jet'};
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
% resampling_method_names_plot_short = {'R-P', 'R-R', 'A-R'};
resampling_method_names_plot_short = {'A-R'};
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
max_error_threshold = 1; %00;
num_bins = 30;
max_error_plot = 0.2;
% bins = linspace(-max_error_threshold, max_error_threshold, num_bins*2);
% bins_abs = linspace(1e-3, max_error_threshold, num_bins);
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
displacement_color_min = [0, 0, 0, -2, 0];
displacement_color_max = [15, 5, 5, 2, 5]; 

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

% ============================
%% loop through window resolutions
% ============================
for window_resolution_index = 1:num_window_resolution
    fprintf('WS: %d\n', window_resolution_index);

    % h1 = figure(1);
    % h2 = figure(2);
    
    % ============================
    %% loop through datasets
    % ============================
    for dataset_index = [3, 4] %1:num_datasets
        % dataset name
        dataset_name = dataset_name_array{dataset_index};
        fprintf('dataset: %s\n', dataset_name);

        % ================================================
        % load results
        % ================================================
        % directory to save results for this case
        current_read_directory = fullfile(top_write_directory, dataset_name, ['WS' num2str(window_resolution_index)], 'weights-rms-change-study-profile');
        % create directory to save figures
        if save_figures
            figure_save_directory = fullfile(current_read_directory, ['figures-' metric_name '-thresh=' num2str(max_error_threshold, '%.2f') 'pix']);
            mkdir_c(figure_save_directory);
        end

        % --------------------------
        % load processing results
        % --------------------------
        % file name
        filename = ['resampling-statistics-nr='  num2str(num_resampling_trials, '%d') '.mat'];
        % load results to file
        load(fullfile(current_read_directory, filename));
        % --------------------------
        % load weights, errors and uncertainties
        % --------------------------
        % file name
        filename = ['resampling-statistics-nr='  num2str(num_resampling_trials, '%d') '-weights-unc-' metric_name '-03.mat'];
        % load results to file
        load(fullfile(current_read_directory, filename));
        % --------------------------
        % load statistics
        % --------------------------
        % file name
        filename = ['resampling-statistics-nr='  num2str(num_resampling_trials, '%d') '-statistics-' metric_name '-new-thresh=' num2str(max_error_threshold, '%.2f') 'pix.mat'];
        % load results to file
        load(fullfile(current_read_directory, filename));

        num_snapshots = size(err.x, 1);
        num_grid_points_prof = size(err.x, 2);
        num_points = numel(err.x);

        % ==========================
        % plot contour and profile line
        % ==========================
        figure
        display_prana_displacements_contour(results_all{window_resolution_index, dataset_index}{1}, [0, displacement_color_max(dataset_index)], 'parula')
        hold on
        % extract co-ordinates
        x = results_all{window_resolution_index, dataset_index}{1}.X(r_prof, c_prof);
        y = results_all{window_resolution_index, dataset_index}{1}.Y(r_prof, c_prof);
        % plot slice line
        plot(x(1, :), y(:, 1), '--', 'color', colors(2, :))        
        title(dataset_name_array_plot{dataset_index})
        set(gcf, 'resize', 'off');
        set(gcf, 'units', 'inches', 'position', [440   355   560   420]/user_screen_resolution);
        set(gcf, 'resize', 'off');
        drawnow();

        if save_figures
            save_figure_to_png_svg_fig(figure_save_directory, 'flow-field', [1, 0, 0])
            % close all
        end

        continue;
        % ==========================
        % plot spatial variation of rms
        % ==========================
        figure        
        for component_index = 1:num_components
            % subplot_index = (dataset_index - 1) * num_components + component_index;

            % component name
            component_name = components{component_index};
            % subplot(num_datasets, num_components, subplot_index)
            subplot(1, num_components, component_index)
            % x axes
            x = linspace(0, 1, num_grid_points_prof);
            % rms error
            % plot(x, err_rms.(component_name), 's', 'color', [0, 0, 0]);
            plot(x, err_rms.(component_name), 'color', [0, 0, 0]);
            hold on
            for method_index = 1:num_individual_methods
                plot(x, unc_indiv_rms.(component_name)(:, method_index), line_symbols{method_index}, 'color', colors(1, :));
                % plot(x, unc_indiv_rms.(component_name)(:, method_index), symbols{method_index}, 'color', colors(1, :));
            end        
            plot(x, unc_comb_rms.(component_name)(:, resampling_method_index_plot), 'color', colors(2, :));
            % plot(x, unc_comb_rms.(component_name)(:, resampling_method_index_plot), symbols{resampling_method_index_plot}, 'color', colors(2, :));
            ylim([0 0.2])

            if dataset_index == 1
                title(upper(component_name))
            end
            box off

            if dataset_index == 3 && component_index == 1
                ylabel('RMS Error and Uncertainty (pix).')
            end
        end

        set_subplots_height(gcf, 0.7);
        
        lgd = legend({'Error', individual_method_array{:}, resampling_method_names_plot_short{:}}, ...
                        'location', 'northoutside', 'orientation', 'horizontal');
        set(gcf, 'resize', 'off');
        set(gcf, 'units', 'inches', 'position', [260 371 787 361]/user_screen_resolution);
        set(gcf, 'resize', 'off');
        drawnow();
        lgd.Position(1:2) = [0.3, 0.9];

        if save_figures
            save_figure_to_png_svg_fig(figure_save_directory, 'rms-err-unc', [1, 0, 0])
        end

        % ==========================
        % plot spatial variation of weights
        % ==========================
        figure
        for component_index = 1:num_components
            % subplot_index = (dataset_index - 1) * num_components + component_index;
            % component name
            component_name = components{component_index};
            % subplot(num_datasets, num_components, subplot_index)
            subplot(1, num_components, component_index)
            % x axes
            x = linspace(0, 1, num_grid_points_prof);
            hold on
            for method_index = 1:num_individual_methods                
                plot(x, wt_med{resampling_method_index_plot}.(component_name)(:, method_index), line_symbols{method_index}, 'color', colors(1, :));
            end        
            ylim([0 1])
            if dataset_index == 1
                title(upper(component_name))
            end
            box off

            if dataset_index == 3 && component_index == 1
                ylabel('Median Weights')
            end
        end

        set_subplots_height(gcf, 0.7);
        
        lgd = legend({individual_method_array{:}}, ...
                        'location', 'northoutside', 'orientation', 'horizontal');
        set(gcf, 'resize', 'off');
        set(gcf, 'units', 'inches', 'position', [260 371 787 361]/user_screen_resolution);
        set(gcf, 'resize', 'off');
        drawnow();
        lgd.Position(1:2) = [0.4, 0.9];        

        if save_figures
            save_figure_to_png_svg_fig(figure_save_directory, 'weights', [1, 0, 0])
        end        
    end    
end
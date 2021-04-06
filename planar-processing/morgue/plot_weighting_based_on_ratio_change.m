% This script performs resampling calculations for a range of particle removal
% percentage and plots the change in uncertainties for different methods

%%
clear
close all
clc

restoredefaultpath;
addpath(genpath('/scratch/shannon/c/aether/Projects/BOS/general-codes/matlab-codes/'));
addpath(genpath('/scratch/shannon/c/aether/Projects/BOS/general-codes/piv-image-generation/'));
% addpath('../prana-uncertainty-average-dc-new-im-2-save-deform-cs/');
addpath ../prana/
addpath ../general-codes/
addpath ../histogram_distance/
setup_default_settings;

% ============================
%% read/write settings
% ============================
% window resolution
window_resolution_array = [64]; %, 32];
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
intensity_threshold = 10;
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

% displacement_contour_levels = linspace(displacement_color_min, displacement_color_max, 100);

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

% start timer
tic
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
        % create directory to save figures
        if save_figures
            figure_save_directory = fullfile(current_read_directory, ['figures-' metric_name]);
            mkdir_c(figure_save_directory);
        end

        % file name
        filename = ['resampling-statistics-nt=' num2str(num_trials, '%d') 'nr='  num2str(num_resampling_trials, '%d') '.mat'];
        % load results to file
        load(fullfile(current_read_directory, filename));
        % file name
        filename = ['resampling-statistics-nt=' num2str(num_trials, '%d') 'nr='  num2str(num_resampling_trials, '%d') '-weights-unc-' metric_name '.mat'];
        % load results to file
        load(fullfile(current_read_directory, filename));

        % directory to save results for this case
        current_read_directory = fullfile(top_write_directory, dataset_name, ['WS' num2str(window_resolution_index)], 'weights-rms-change-study');
        % create directory to save figures
        if save_figures
            figure_save_directory = fullfile(current_read_directory, ['figures-' metric_name]);
            mkdir_c(figure_save_directory);
        end
        save_figures = 1;

        % % ==========================
        % % plot variation of  with particle % for just one trial
        % % ==========================
        % trial_index = 100;
        % ppr = percentage_particles_remove_array*100;

        % plot_unc_ratio_with_particle_removal(unc_ratio, unc_ratio_fit, trial_index, ppr, ...
        %                                     individual_method_array, resampling_method_names_plot, ['Uncertainty Ratio ' upper(metric_name)], user_screen_resolution);

        % if save_figures
        %     save_figure_to_png_svg_fig(figure_save_directory, 'rate', [1, 0, 0])
        % end

        % % % ==========================
        % % % compare error models
        % % % ==========================
        % % compare_error_models(error_models, err_trials, err_est_indiv, err_est_comb, individual_method_array, resampling_method_names_plot_short, ...
        % %                         bins, user_screen_resolution);
        % % if save_figures
        % %     save_figure_to_png_svg_fig(figure_save_directory, 'error-model-comparison', [1, 0, 0]);
        % % end

        % % ==========================
        % % plot pdf of weights
        % % ==========================
        % plot_pdf_weights(bins_w, pdf_w, individual_method_array, resampling_method_names_plot_short, colors, line_symbols, user_screen_resolution)

        % if save_figures
        %     save_figure_to_png_svg_fig(figure_save_directory, 'pdf-weights', [1, 0, 0]);
        % end

        % ==========================
        % plot pdfs of errors and uncertainties
        % ==========================
        unc_indiv_all = [unc_indiv.x; unc_indiv.y];
        unc_comb_all = [unc_comb.x; unc_comb.y];
        unc_indiv_2 = mat2cell(unc_indiv_all, num_trials*2, ones(1, num_individual_methods));
        unc_comb_2 = mat2cell(unc_comb_all, num_trials*2, ones(1, num_resampling_methods));
        sigma_rms_indv_2 = nanrms(unc_indiv_all, 1);
        sigma_rms_comb_2 = nanrms(unc_comb_all, 1);
        violins = make_violin_plot_error_uncertainty_02(err_trials, unc_indiv_2, unc_comb_2, ...
                                                        err_rms, sigma_rms_indv_2, sigma_rms_comb_2, ...
                                                        bins, individual_method_array, resampling_method_names_plot_short, ...
                                                        colors, user_screen_resolution, max_error_threshold);
        
        % face colors
        color_all = cell(1, numel(violins));

        % add lines corresponding to rms
        violin_index = 4;
        % x co-ordinates of current violin plot
        x = violins(violin_index).XData;
        % plot confidence intervals for the error
        plot([min(violins(1).XData), max(violins(numel(violins)).XData)], confint_err*[1, 1], '-*', 'color', [0, 0, 0, 0.2])
        % plot([min(violins(1).ViolinPlot.XData), max(violins(numel(violins)).ViolinPlot.XData)], confint_err_abs*[1, 1], '-o', 'color', [0, 0, 0, 0.2])

        if save_figures
            save_figure_to_png_svg_fig(figure_save_directory, 'violin-mean-sub-confint-new', [1, 0, 0]);
        end

        % ===================================
        %% comparison of true vs estimated error
        % ===================================
        for model_index = 1:num_models
            error_model = error_models{model_index};
            % ===================================
            %% quantile-quantile plot 
            % ===================================
            figure
            qq_plot_err_true_est(err_trials, {err_est_indiv{:, model_index}}, {err_est_comb{:, model_index}}, ...
                                individual_method_array, resampling_method_names_plot_short, ...
                                max_error_threshold, 2, colors, symbols, user_screen_resolution);

            drawnow();

            % save figures
            if save_figures
                save_figure_to_png_svg_fig(figure_save_directory, ['qq-error-true-vs-est-' error_model], [1, 0, 0]);
                % save_figure_to_png_svg_fig(figure_save_directory, ['qq-error-true-vs-est-' error_pdf], [1, 0, 0]);
                % export_fig(fullfile(figure_write_directory, 'qq-error-true-vs-est.png'), '-r600');
            end

            % ===================================
            % total variation distance 
            % ===================================
            figure
            % create array of distances
            Y = [d_err_est_individual(:, model_index)', d_err_est_combined(:, model_index)'];
            make_histogram_plot(histogram_distance_categories, Y, colors);
            drawnow();
            % save figure
            if save_figures
                save_figure_to_png_svg_fig(figure_save_directory, ['histogram-distance-' error_model], [1, 0, 0]);
                % save_figure_to_png_svg_fig(figure_save_directory, ['histogram-distance-' error_pdf '-abs'], [1, 0, 0]);
            else
                pause(0.1);
            end    
        end
    end
end

% stop timer
toc

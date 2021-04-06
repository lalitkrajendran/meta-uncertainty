clear
close all
clc

restoredefaultpath;
addpath(genpath('/scratch/shannon/c/aether/Projects/BOS/general-codes/matlab-codes/'));
setup_default_settings;
% addpath('../prana-uncertainty-average-dc-new-im-2-save-deform-cs/');
addpath('../prana/');
addpath('../general-codes/')

% ========================
%% read/write settings
% ========================
% window resolution
window_resolution_array = [64, 32];
num_window_resolution = numel(window_resolution_array);

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
% combination_method_array = {'unwt'; 'var-covar'; 'entropy'}; % 'prob'};
% combination_method_array = {'unwt'; 'pd-var'; 'entropy'}; % 'prob'};
combination_method_array = {'unwt'; 'var'; 'entropy'}; % 'prob'};
num_combination_methods = numel(combination_method_array);
% combination_method_names = {'Unwt'; 'Var-Covar'; 'Entropy'};
% combination_method_names = {'Unwt'; 'PD-Var'; 'Entropy'};
combination_method_names = {'Unwt'; 'Var'; 'Entropy'};

% methods to calculate histogram distances
histogram_distance_methods = {'total_variation_distance'; 'chi_square_statistics'; 'kolmogorov_smirnov_distance'; ...
    'hellinger_distance'}; % 'kullback_leibler_divergence'};
num_distance_methods = numel(histogram_distance_methods);

histogram_distance_categories = categorical({individual_method_array{:}, combination_method_names{:}});
% sort in the desired order
histogram_distance_categories = reordercats(histogram_distance_categories, ...
                                {individual_method_array{:}, combination_method_names{:}});

% directory where results of this analysis are to be saved
top_write_directory = fullfile('/scratch/shannon/a/aether/Projects/PIV_PTV_Uncertainty_Quantification/uncertainty-combination/analysis/results/planar-dataset', 'monte-carlo', 'individual-datasets');
mkdir_c(top_write_directory);

% ========================
%% statistical analysis settings
% ========================
% number of trials
num_trials = 1e3;
% minimum allowable error (pix.)
min_error_threshold = 1e-3;
% maximum allowable error (pix.)
max_error_threshold = 0.2;
% number of bins for histogram
% num_bins = 30;
num_bins = round(max_error_threshold/min_error_threshold * 0.2);
% bins for histograms
% bins = linspace(min_error_threshold, max_error_threshold, num_bins);
bins = linspace(min_error_threshold, max_error_threshold, num_bins);
% number of bins for the coarse histogram
num_bins_coarse = 8;
% bins for histogram of weights
bins_weights = linspace(0, 1, 25);

% ========================
%% resampling settings
% ========================
% number of particles to remove
% num_particles_remove_array = [5, 10, 15];
percentage_particles_remove = 0.25;
% number of resampling trials
num_resampling_trials = 1e3;
% uncertainty flags for the various methods (0 = False, 1 = True)
uncertainty_flags.ppruncertainty = 0;
uncertainty_flags.miuncertainty = 0;
uncertainty_flags.mcuncertainty = 1;

resampling_case_name = ['ppr' num2str(percentage_particles_remove, '%.2f') '_nrs' num2str(num_resampling_trials) '-prewindow-randpart-snr-newproc'];
% ========================
%% plot settings
% ========================
% user screen resolution
user_screen_resolution = 113;
% save_figure? (true/false)
save_figures = 1;
% symbols
symbols = {'o', '^', 'v'};
% colors
colors = lines(3);
% line symbols
line_symbols = {'-'; ':'; '-.'};
% number of data points to skip for the qq plot
num_skip = 10;
% number of trials to skip for the pdf plot
skip_trials = 100;

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
    
        % ===============================
        %% directory settings for this case
        % ===============================
        % directory to load results for this case
        current_read_directory = fullfile(top_write_directory, dataset_name_array{dataset_index}, ['WS' num2str(window_resolution_index)], ...
        ['trials=' num2str(num_trials, '%d')], resampling_case_name, ...
        ['max_error=' num2str(max_error_threshold, '%.2f') 'pix']);

        figure_write_directory = fullfile(current_read_directory, 'figures-var');
        mkdir_c(figure_write_directory);

        %% load results
        % load(fullfile(current_read_directory, 'error_uncertainty_statistics.mat'));
        load(fullfile(current_read_directory, 'error_uncertainty_statistics_02.mat'));

        return;
        % =============================
        %% plot unc sub vs unc
        % =============================
        figure
        for method_index = 1:num_individual_methods
            subplot(1, num_individual_methods, method_index)

            plot(unc_individual_all{method_index}, unc_individual_sub_all{method_index}, 'o')
            hold on
            plot([0 1], [0 1], 'k')
            box off
            axis equal
            axis([0 0.2 0 0.2])
            xlabel('Uncertainty, Original (pix.)')
            ylabel('Uncertainty, Sub (pix.)')
            title(individual_method_array{method_index})
        end

        set(gcf, 'resize', 'off')
        set(gcf, 'position', [202   378   840   400])
        drawnow();
        % save figure
        if save_figures
            save_figure_to_png_svg_fig(figure_write_directory, 'unc-sub-vs-orig', [1, 0, 0]);
            % export_fig(fullfile(figure_write_directory, 'unc-vs-snr-mean.png'), '-r600'); 
        else
            pause(0.1);
        end
        return;
        % =============================
        %% plot snr metric for resampled uncertainties
        % =============================
        figure
        for snr_method_index = 1:num_snr_methods
            subplot(1, num_snr_methods, snr_method_index)
            for uncertainty_method_index = 1:num_individual_methods
                plot(mean_snr_resampled(snr_method_index, :), mean_unc_resampled(uncertainty_method_index, :), 'o')
                hold on
            end

            box off
            xlabel(snr_metric_array{snr_method_index})
            ylabel('Uncertainty (pix.)')
            legend(individual_method_array, 'location', 'northoutside', 'orientation', 'horizontal')
            % legend boxon
        end

        set(gcf, 'resize', 'off')
        set(gcf, 'position', [202   378   840   400])
        drawnow();
        % save figure
        if save_figures
            save_figure_to_png_svg_fig(figure_write_directory, 'error-uncertainty-histograms-violin', [1, 0, 0]);
            % export_fig(fullfile(figure_write_directory, 'unc-vs-snr-mean.png'), '-r600'); 
        else
            pause(0.1);
        end

        % =============================
        %% violin plot of error and uncertainty histograms
        % =============================
        make_violin_plot_error_uncertainty(abs(err_all_valid'), unc_individual_all_valid, unc_combined_all_valid, rms_err, rms_unc_individual, rms_unc_combined, ...
                                            individual_method_array, combination_method_array, ...
                                            colors, user_screen_resolution, max_error_threshold);

        drawnow();
        % save figure
        if save_figures
            save_figure_to_png_svg_fig(figure_write_directory, 'error-uncertainty-histograms-violin', [1, 0, 0]);
            % export_fig(fullfile(figure_write_directory, 'error-uncertainty-histograms-violin.png'), '-r600'); 
        else
            pause(0.1);
        end

        % =============================
        %% individual vs resampled uncertaintes
        % =============================
        plot_pdf_uncertainties_individual_resampled(pdf_unc_individual, pdf_unc_resampled, ...
                                                    rms_unc_individual, rms_unc_resampled, ...
                                                    bins, individual_method_array, user_screen_resolution);
        drawnow();
        % save figure
        if save_figures
            save_figure_to_png_svg_fig(figure_write_directory, 'pdf-uncertainty-individual-resampled', [1, 0, 0]);
            % export_fig(fullfile(figure_write_directory, 'error-uncertainty-histograms-violin.png'), '-r600'); 
        else
            pause(0.1);
        end
        
        % =============================
        %% resampled uncertainties for each trial
        % =============================
        plot_pdf_error_uncertainties_resampled_trials(pdf_err_resampled_trials, pdf_unc_resampled_trials, bins, individual_method_array, skip_trials, user_screen_resolution)
        drawnow();

        % save figure
        if save_figures
            save_figure_to_png_svg_fig(figure_write_directory, ['pdf-uncertainty-resampled-trials-skip=' num2str(skip_trials)], [1, 0, 0]);
            % export_fig(fullfile(figure_write_directory, 'error-uncertainty-histograms-violin.png'), '-r600'); 
        else
            pause(0.1);
        end

        % ========================
        %% pdf of weights
        % ========================
        plot_pdf_weights(bins_weights, pdf_w, individual_method_array, combination_method_array, colors, line_symbols, user_screen_resolution);
        drawnow();        
        % save figures
        if save_figures
            save_figure_to_png_svg_fig(figure_write_directory, 'pdf-weights', [1, 0, 0]);
            % export_fig(fullfile(figure_write_directory, 'pdf-weights.png'), '-r600');
        end

        % ===================================
        %% quantile-quantile plot of true vs estimated error
        % ===================================
        figure
        qq_plot_err_true_est(err_all_valid, err_est_valid_individual, err_est_valid_combined, ...
                            individual_method_array, combination_method_array, max_error_threshold, num_skip, colors, symbols, user_screen_resolution);
        drawnow();
        % save figures
        if save_figures
            save_figure_to_png_svg_fig(figure_write_directory, 'qq-error-true-vs-est', [1, 0, 0]);
            % export_fig(fullfile(figure_write_directory, 'qq-error-true-vs-est.png'), '-r600');
        end

        % ===================================
        % plot pdf of true vs estimated error
        % ===================================
        plot_pdf_err_true_est(bins, pdf_err, pdf_err_est_individual, pdf_err_est_combined, max_error_threshold, individual_method_array, combination_method_names, ...
                            colors, line_symbols, user_screen_resolution);

        drawnow();
        % save figures        
        if save_figures
            save_figure_to_png_svg_fig(figure_write_directory, 'pdf-error-true-vs-est', [1, 0, 0]);
            % export_fig(fullfile(figure_write_directory, 'pdf-error-true-vs-est.png'), '-r600');
        end
    
        % ===================================
        % plot cdf of true vs estimated error
        % ===================================
        plot_cdf_err_true_est(bins, cdf_err, cdf_err_est_individual, cdf_err_est_combined, max_error_threshold, individual_method_array, combination_method_names, ...
                            colors, line_symbols, user_screen_resolution);
        
        drawnow();
        % save figures
        if save_figures
            save_figure_to_png_svg_fig(figure_write_directory, 'cdf-error-true-vs-est', [1, 0, 0]);
            % export_fig(fullfile(figure_write_directory, 'cdf-error-true-vs-est.png'), '-r600');
        end

        % ===================================
        % plot total variation distance
        % ===================================
        figure
        distance_method_index = 1;
        % create array of distances
        Y = [d_err_est_individual(:, distance_method_index)', d_err_est_combined(:, distance_method_index)'];
        make_histogram_plot(histogram_distance_categories, Y, colors);
        drawnow();
        % save figure
        if save_figures
            save_figure_to_png_svg_fig(figure_write_directory, 'histogram-distance', [1, 0, 0]);
        else
            pause(0.1);
        end
        return;
    end
end

% =============================
%% make plots for all datasets combined
% =============================
% directory to load results for this case
current_read_directory = fullfile(top_write_directory, 'all-combined', ...
['trials=' num2str(num_trials, '%d')], resampling_case_name, ...
['max_error=' num2str(max_error_threshold, '%.2f') 'pix-new']);

% load aggregated results
% load(fullfile(current_read_directory, 'error_uncertainty_statistics.mat'));
load(fullfile(current_read_directory, 'error_uncertainty_statistics_02.mat'));

figure_write_directory = fullfile(current_read_directory, 'figures-var');
mkdir_c(figure_write_directory);

% =============================
%% violin plot of error and uncertainty histograms
% =============================
% make_violin_plot_error_uncertainty(abs(err_all_valid)', unc_individual_all_valid{1}', unc_individual_all_valid{2}', unc_individual_all_valid{3}', ...
% unc_combined_all_valid{1}', unc_combined_all_valid{2}', unc_combined_all_valid{3}');
make_violin_plot_error_uncertainty(abs(err_all_valid_all_datasets), sigma_all_valid_all_datasets, unc_combined_all_valid_all_datasets, ...
                                        err_rms_all_datasets, sigma_rms_all_datasets, sigma_rms_comb_all_datasets, ...
                                        individual_method_array, combination_method_array, ...
                                        colors, user_screen_resolution, max_error_threshold);

drawnow();
% save figure
if save_figures
    % save_figure_to_png_svg_fig(figure_write_directory, 'error-uncertainty-histograms-violin', [1, 0, 0]);
    export_fig(fullfile(figure_write_directory, 'error-uncertainty-histograms-violin.png'), '-r600'); 
else
    pause(0.1);
end


% =============================
%% individual vs resampled uncertaintes
% =============================
plot_pdf_uncertainties_individual_resampled(pdf_unc_individual_all_datasets, pdf_unc_resampled_all_datasets, ...
                                            rms_unc_individual_all_datasets, rms_unc_resampled_all_datasets, ...
                                            bins, individual_method_array, user_screen_resolution);

drawnow();
% save figure
if save_figures
    save_figure_to_png_svg_fig(figure_write_directory, 'pdf-uncertainty-individual-resampled', [1, 0, 0]);
    % export_fig(fullfile(figure_write_directory, 'error-uncertainty-histograms-violin.png'), '-r600'); 
else
    pause(0.1);
end

% ===================================
%% quantile-quantile plot of true vs estimated error
% ===================================
figure
qq_plot_err_true_est(err_all_valid_all_datasets, err_est_valid_individual_all_datasets, err_est_valid_combined_all_datasets, ...
                    individual_method_array, combination_method_array, max_error_threshold, num_skip * 10, colors, symbols, user_screen_resolution);
drawnow();
% save figures
if save_figures
    % save_figure_to_png_svg_fig(figure_write_directory, 'qq-error-true-vs-est', [1, 0, 0]);
    export_fig(fullfile(figure_write_directory, 'qq-error-true-vs-est.png'), '-r600');
end

% ===================================
% plot total variation distance
% ===================================
figure
distance_method_index = 1;
% create array of distances
Y = [d_err_est_individual_all_datasets(:, distance_method_index)', d_err_est_combined_all_datasets(:, distance_method_index)'];
make_histogram_plot(histogram_distance_categories, Y, colors);
drawnow();
% save figure
if save_figures
    save_figure_to_png_svg_fig(figure_write_directory, 'histogram-distance', [1, 0, 0]);
else
    pause(0.1);
end

% ========================
%% pdf of weights
% ========================
plot_pdf_weights(bins_weights, pdf_w_all_datasets, individual_method_array, combination_method_array, colors, line_symbols, user_screen_resolution);
drawnow();        
% save figures
if save_figures
    % save_figure_to_png_svg_fig(figure_write_directory, 'pdf-weights', [1, 0, 0]);
    export_fig(fullfile(figure_write_directory, 'pdf-weights.png'), '-r600');
end

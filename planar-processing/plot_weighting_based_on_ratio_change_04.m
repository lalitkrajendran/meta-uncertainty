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
dataset_name_array = {'PivChal03B'}; % 'PivChal05B'; 'stagnation_flow'; 'Vortex_Ring'; 'Jetdata'};
% dataset_name_array = {'Vortex_Ring'}; 
dataset_name_array_plot = {'TBL'; 'LSB'; 'SF'; 'VR'; 'Jet'};
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
ppr = percentage_particles_remove_array*100;
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
displacement_color_min = [0, 0, 0, 0, 0];
displacement_color_max = [15, 5, 5, 5, 5]; 

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

% % ============================
% %% loop through window resolutions
% % ============================
% for window_resolution_index = 1:num_window_resolution
%     fprintf('WS: %d\n', window_resolution_index);

%     % ============================
%     %% loop through datasets
%     % ============================
%     for dataset_index = 1:num_datasets
%         % dataset name
%         dataset_name = dataset_name_array{dataset_index};
%         fprintf('dataset: %s\n', dataset_name);

%         % ================================================
%         % load results
%         % ================================================
%         % directory to save results for this case
%         current_read_directory = fullfile(top_write_directory, dataset_name, ['WS' num2str(window_resolution_index)], 'weights-rms-change-study');
%         % create directory to save figures
%         if save_figures
%             figure_save_directory = fullfile(current_read_directory, ['figures-' metric_name '-thresh=' num2str(max_error_threshold, '%.2f') 'pix']);
%             mkdir_c(figure_save_directory);
%         end

%         % load errors
%         errors = errors_all{window_resolution_index, dataset_index};

%         % --------------------------
%         % load processing results
%         % --------------------------
%         % file name
%         filename = ['resampling-statistics-nt=' num2str(num_trials, '%d') 'nr='  num2str(num_resampling_trials, '%d') '.mat'];
%         % load results to file
%         load(fullfile(current_read_directory, filename));
%         % --------------------------
%         % load weights, errors and uncertainties
%         % --------------------------
%         % file name
%         filename = ['resampling-statistics-nt=' num2str(num_trials, '%d') 'nr='  num2str(num_resampling_trials, '%d') '-weights-unc-' metric_name '-02-new.mat'];
%         % load results to file
%         load(fullfile(current_read_directory, filename));
%         % --------------------------
%         % load statistics
%         % --------------------------
%         % file name
%         % filename = ['resampling-statistics-nt=' num2str(num_trials, '%d') 'nr='  num2str(num_resampling_trials, '%d') '-statistics-' metric_name '-new.mat'];
%         filename = ['resampling-statistics-nt=' num2str(num_trials, '%d') 'nr='  num2str(num_resampling_trials, '%d') '-statistics-' metric_name '-new-thresh=' num2str(max_error_threshold, '%.2f') 'pix.mat'];
%         % load results to file
%         load(fullfile(current_read_directory, filename));
%         % --------------------------
%         % load physical properties
%         % --------------------------
%         % file name
%         filename = ['resampling-statistics-nt=' num2str(num_trials, '%d') 'nr='  num2str(num_resampling_trials, '%d') '-physical-prop.mat'];
%         % load results
%         load(fullfile(current_read_directory, filename));

%         % ================================================
%         %% determine indices to be used for this dataset
%         % ================================================
%         if strcmpi(dataset_name, 'Vortex_Ring')        
%             % extract region with non-nan values
%             mask = ~isnan(errors.err_U(:, :, 1));
%             [cmin, cmax, rmin, rmax] = extract_image_mask_extents(mask(:, :, 1));
            
%             % % extract co-ordinate arrays for the processing grid
%             % xvec = errors.X(1, :);
%             % yvec = errors.Y(:, 1);
            
%             % find vector locations in the processing results for which
%             % true solution is available
%             r = rmin:rmax;
%             c = cmin:cmax;
%             % c = find(xvec > errors.X(1, 1) & xvec < errors.X(1, end));
%             % r = find(yvec > errors.Y(1, 1) & yvec < errors.Y(end, 1));
%         elseif strcmp(dataset_name, 'Jetdata')
%             % index limits for processing
%             % X 320 to 392 (72); Y 205 to 325 (120)
%             cmin_p = 5;
%             cmax_p = 23; %303+18-1=320  303+90-1=392  X direction
%             rmin_p = 5;
%             rmax_p = 35; %188+18-1=205  188+138-1=325 Y direction
                        
%             % only retain true solution in the index limits for
%             % processing
%             c = find(errors.X(1, :) >= 320  & errors.X(1, :) <= 392);
%             r = find(errors.Y(:, 1) >= 205  & errors.Y(:, 1) <= 325);
%         else
%             c = 1:size(errors.X, 2);
%             r = 1:size(errors.X, 1);
%         end

%         % --------------------------
%         % extract boundaries
%         % --------------------------
%         xmin = results_all{window_resolution_index, dataset_index}{1}.X(1, c(1));
%         xmax = results_all{window_resolution_index, dataset_index}{1}.X(1, c(end));
%         ymin = results_all{window_resolution_index, dataset_index}{1}.Y(r(1), 1);
%         ymax = results_all{window_resolution_index, dataset_index}{1}.Y(r(end), 1);

%         % % ==========================
%         % % plot contour
%         % % ==========================
%         % figure
%         % display_prana_displacements_contour(results_all{window_resolution_index, dataset_index}{1}, [0, displacement_color_max(dataset_index)], 'jet')
%         % axis([xmin, xmax, ymin, ymax])
%         % title(dataset_name_array_plot{dataset_index})
%         % set(gcf, 'resize', 'off');
%         % set(gcf, 'units', 'inches', 'position', [440   355   560   420]/user_screen_resolution);
%         % set(gcf, 'resize', 'off');
%         % drawnow();

%         % if save_figures
%         %     save_figure_to_png_svg_fig(figure_save_directory, 'flow-field', [1, 0, 0])
%         %     % close all
%         % end

%         % % ==========================
%         % % plot uncertainty as function of fractional displacement
%         % % ==========================
%         % frac_disp_U = U_trials_valid - round(U_trials_valid);
%         % frac_disp_V = V_trials_valid - round(V_trials_valid);
%         % x = [frac_disp_U; frac_disp_V];

%         % frac_disp_U = U_ref_trials_valid - round(U_ref_trials_valid);
%         % frac_disp_V = V_ref_trials_valid - round(V_ref_trials_valid);
%         % x_ref = [frac_disp_U; frac_disp_V];

%         % num_bins_disp = 10;
%         % bins_disp = linspace(-0.5, 0.5, num_bins_disp);
%         % Y = discretize(x, bins_disp);
%         % Y_ref = discretize(x_ref, bins_disp);
%         % err = [err_rms_U_valid; err_rms_V_valid];

%         % err_rms_binned = nans(num_bins_disp, 1);
%         % err_rms_binned_ref = nans(num_bins_disp, 1);
        
%         % unc_indiv_rms_binned = nans(num_bins_disp, num_individual_methods);
%         % unc_indiv_rms_binned_ref = nans(num_bins_disp, num_individual_methods);

%         % unc_comb_rms_binned = nans(num_bins_disp, num_individual_methods);
%         % unc_comb_rms_binned_ref = nans(num_bins_disp, num_individual_methods);

%         % for bin_index = 1:num_bins_disp
%         %     indices = find(Y == bin_index);
%         %     indices_ref = find(Y_ref == bin_index);

%         %     err_rms_binned(bin_index) = nanrms(err(indices)');
%         %     err_rms_binned_ref(bin_index) = nanrms(err(indices_ref)');

%         %     unc_indiv_rms_binned(bin_index, :) = nanrms(unc_indiv_all(indices, :), 1);
%         %     unc_indiv_rms_binned_ref(bin_index, :) = nanrms(unc_indiv_all(indices_ref, :), 1);

%         %     unc_comb_rms_binned(bin_index, :) = nanrms(unc_comb_all(indices, :), 1);
%         %     unc_comb_rms_binned_ref(bin_index, :) = nanrms(unc_comb_all(indices_ref, :), 1);
%         % end

%         % figure
%         % subplot(1, 2, 1)
%         % plot(bins_disp, err_rms_binned, 'ks')
%         % hold on
%         % for method_index = 1:num_individual_methods
%         %     plot(bins_disp, unc_indiv_rms_binned(method_index), symbols{method_index}, 'color', colors(1, :))
%         % end
%         % plot(bins_disp, unc_comb_rms_binned(resampling_method_index_plot), symbols{resampling_method_index_plot}, 'color', colors(2, :))
%         % box off
%         % ylim([0 0.2])
%         % xlabel('Fractional Displacement (pix.)')
%         % title('Measured Displacement')

%         % subplot(1, 2, 2)
%         % plot(bins_disp, err_rms_binned_ref, 'ks')
%         % hold on
%         % for method_index = 1:num_individual_methods
%         %     plot(bins_disp, unc_indiv_rms_binned_ref(method_index), symbols{method_index}, 'color', colors(1, :))
%         % end
%         % plot(bins_disp, unc_comb_rms_binned_ref(resampling_method_index_plot), symbols{resampling_method_index_plot}, 'color', colors(2, :))
%         % box off
%         % ylim([0 0.2])
%         % xlabel('Fractional Displacement (pix.)')
%         % title('True Displacement')        
%         % lgd = legend({'Error', individual_method_array{:}, resampling_method_names_plot_short{:}}, ...
%         %                     'location', 'northoutside', 'orientation', 'horizontal');
%         % set_subplots_height(gcf, 0.7);

%         % set(gcf, 'resize', 'off')
%         % set(gcf, 'units', 'inches', 'position', [112         348        1138         735]/user_screen_resolution);
%         % set(gcf, 'resize', 'off')
%         % drawnow();
%         % lgd.Position(1:2) = [0.35, 0.9];
%         % return;

%         % return;
%         % % ==========================
%         % % plot pdf of unc ratio with particle % 
%         % % ==========================
%         % trial_index = 100;
%         % plot_unc_ratio_violin({unc_resampling_trials{trial_index, :, resampling_method_index_plot}}, unc_sub_trials{trial_index}, ppr, num_resampling_trials, ...
%         %                         individual_method_array, user_screen_resolution)
%         % if save_figures
%         %     save_figure_to_png_svg_fig(figure_save_directory, 'unc-ratio-pdf', [1, 0, 0])
%         % end

%         % return;
%         % % ==========================
%         % % plot pdf of snr ratio with particle % 
%         % % ==========================
%         % trial_index = 100;
%         % plot_snr_ratio_violin({snr_resampling_trials{trial_index, :, resampling_method_index_plot}}, snr_sub_trials{trial_index}, ppr, num_resampling_trials, ...
%         %                         snr_metric_array, user_screen_resolution)

%         % if save_figures
%         %     save_figure_to_png_svg_fig(figure_save_directory, 'snr-ratio-pdf', [1, 0, 0])
%         % end
        
%         % return;
%         % % ==========================
%         % % plot variation of uncertainty with particle % for just one trial
%         % % ==========================
%         % trial_index = 100;
%         % % plot_unc_ratio_with_particle_removal(unc_ratio, unc_ratio_fit, trial_index, ppr, ...
%         % %                                     individual_method_array, resampling_method_names_plot, ['Uncertainty Ratio ' upper(metric_name)], user_screen_resolution);
%         % plot_unc_ratio_with_particle_removal_single(unc_ratio{resampling_method_index_plot}, unc_ratio_fit{resampling_method_index_plot}, wt_trials{resampling_method_index_plot}, ...
%         %                                         trial_index, ppr, individual_method_array, resampling_method_names_plot_short, ...
%         %                                     ['Uncertainty Ratio ' upper(metric_name)], symbols, line_symbols, user_screen_resolution);
        
        
%         % if save_figures
%         %     save_figure_to_png_svg_fig(figure_save_directory, 'rate-new-symbols', [1, 0, 0])
%         % end

%         % return;
%         % % ==========================
%         % % plot pdf of weights
%         % % ==========================
        
%         % % plot_pdf_weights(bins_w, pdf_w, individual_method_array, resampling_method_names_plot_short, colors, line_symbols, user_screen_resolution)
%         % plot_pdf_weights(bins_w, pdf_w{resampling_method_index_plot, :}, individual_method_array, resampling_method_names_plot_short, colors, line_symbols, user_screen_resolution)
%         % resampling_method_index_plot
%         % if save_figures
%         %     save_figure_to_png_svg_fig(figure_save_directory, 'pdf-weights', [1, 0, 0]);
%         % end

%         % % ==========================
%         % % plot pdfs of errors and uncertainties
%         % % ==========================
%         % unc_indiv_2 = mat2cell(unc_indiv_all, num_valid_trials*2, ones(1, num_individual_methods));
%         % % unc_comb_2 = mat2cell(unc_comb_all, num_valid_trials*2, ones(1, num_resampling_methods));
%         % unc_comb_2 = mat2cell(unc_comb_all(:, resampling_method_index_plot), num_valid_trials*2, ones(1, 1)); %num_resampling_methods));
%         % violins = make_violin_plot_error_uncertainty_02(err_all, unc_indiv_2, unc_comb_2, ...
%         %                                                 err_rms, unc_indiv_rms, unc_comb_rms, ...
%         %                                                 bins, individual_method_array, resampling_method_names_plot_short, ...
%         %                                                 colors, user_screen_resolution, max_error_threshold);
                
%         % if save_figures
%         %     save_figure_to_png_svg_fig(figure_save_directory, 'violin', [1, 0, 0]);
%         % end

%         % ===================================
%         %% comparison of true vs estimated error
%         % ===================================
%         error_model = 'gaussian'; 
%         % ===================================
%         %% quantile-quantile plot 
%         % ===================================
%         figure
%         qq_plot_err_true_est(err_all, err_est_indiv, err_est_comb(:, resampling_method_index_plot), ...
%                             individual_method_array, resampling_method_names_plot_short, ...
%                             max_error_plot, 2, colors, symbols, user_screen_resolution);

%         drawnow();

%         return;
%         % save figures
%         if save_figures
%             save_figure_to_png_svg_fig(figure_save_directory, ['qq-error-true-vs-est-' error_model], [1, 0, 0]);
%             % save_figure_to_png_svg_fig(figure_save_directory, ['qq-error-true-vs-est-' error_pdf], [1, 0, 0]);
%             % export_fig(fullfile(figure_write_directory, 'qq-error-true-vs-est.png'), '-r600');
%         end

%         % ===================================
%         % total variation distance 
%         % ===================================
%         figure
%         % create array of distances
%         % Y = [d_err_est_indiv, d_err_est_comb];
%         Y = [d_err_est_indiv, d_err_est_comb(resampling_method_index_plot)];
%         make_histogram_plot(histogram_distance_categories, Y, colors);
%         drawnow();
%         % save figure
%         if save_figures
%             save_figure_to_png_svg_fig(figure_save_directory, ['histogram-distance-' error_model], [1, 0, 0]);
%             % save_figure_to_png_svg_fig(figure_save_directory, ['histogram-distance-' error_pdf '-abs'], [1, 0, 0]);
%         else
%             pause(0.1);
%         end    
        

%         % return;
%     end
% end

% return;
% ==========================
% plot results for all cases combined
% ==========================
% directory to save results for this case
write_directory = fullfile(top_write_directory, 'weights-rms-change-study-consolidated');
% file name
filename = ['resampling-statistics-nt=' num2str(num_trials, '%d') 'nr='  num2str(num_resampling_trials, '%d') '-statistics-' metric_name '-new-thresh=' num2str(max_error_threshold, '%.2f') 'pix.mat'];
% load results
load(fullfile(write_directory, filename));
% file name
filename = ['resampling-statistics-nt=' num2str(num_trials, '%d') 'nr='  num2str(num_resampling_trials, '%d') '-physical-prop.mat'];
% load results
load(fullfile(write_directory, filename));

% create directory to save figures
if save_figures
    figure_save_directory = fullfile(write_directory, ['figures-' metric_name '-thresh=' num2str(max_error_threshold, '%.2f') 'pix']);
    mkdir_c(figure_save_directory);
end


% % ==========================
% % plot results as function of flow properties
% % ==========================
% marker_size = 4;

% figure

% % --------------------------
% % fractional displacement
% % --------------------------
% subplot(2, 2, 1)
% plot(bins_disp, err_rms_binned_disp_ref, 'k')
% hold on
% for method_index = 1:num_individual_methods
%     plot(bins_disp, unc_indiv_rms_binned_disp_ref(:, method_index), symbols{method_index}, 'color', colors(1, :), 'markersize', marker_size, 'linewidth', 1.5)
% end
% plot(bins_disp, unc_comb_rms_binned_disp_ref(:, resampling_method_index_plot), symbols{1}, 'color', colors(2, :), 'markersize', marker_size, 'linewidth', 1.5)
% box off
% xlim([-0.5 0.5])
% ylim([0 0.2])
% % xlabel('Fractional Displacement (pix.)')
% set(gca, 'xticklabel', [])
% ylabel({'RMS Error and'; 'Uncertainty (pix.)'}, 'fontsize', 12)

% subplot(2, 2, 3)
% bar(bins_disp, num_indices_disp_ref/sum(num_indices_disp_ref) * 100)
% box off
% xlim([bins_disp(1) bins_disp(end)])
% ylim([0 75])
% % ylim([0 num_valid_trials_consolidated])
% xlabel('Displacement (pix.)', 'fontsize', 14)
% ylabel('Bin Count (%)', 'fontsize', 14)

% % --------------------------
% % shear
% % --------------------------
% subplot(2, 2, 2)
% plot(bins_grad, err_rms_binned_grad_ref, 'k')
% hold on
% for method_index = 1:num_individual_methods
%     plot(bins_grad, unc_indiv_rms_binned_grad_ref(:, method_index), symbols{method_index}, 'color', colors(1, :), 'markersize', marker_size, 'linewidth', 1.5)
% end
% plot(bins_grad, unc_comb_rms_binned_grad_ref(:, resampling_method_index_plot), symbols{1}, 'color', colors(2, :), 'markersize', marker_size, 'linewidth', 1.5)
% box off
% xlim([bins_grad(1) bins_grad(end)])
% ylim([0 0.2])
% set(gca, 'yticklabel', [])
% set(gca, 'xticklabel', [])

% subplot(2, 2, 4)
% % bar(bins_grad, num_indices_ref/(2 * num_valid_trials_consolidated) * 100)
% bar(bins_grad, num_indices_grad_ref/sum(num_indices_grad_ref) * 100)
% box off
% xlim([bins_grad(1) bins_grad(end)])
% ylim([0 75])
% % ylim([0, 2 * num_valid_trials_consolidated])
% xlabel('Shear (pix./pix.)', 'fontsize', 14)
% % ylabel('Bin Count')
% set(gca, 'yticklabel', [])

% subplot(2, 2, 1)
% lgd = legend({'Error', individual_method_array{:}, 'Comb'}, ...
%                     'location', 'northoutside', 'orientation', 'horizontal', 'fontsize', 12);

% % set_subplots_height(gcf, 0.7);

% set(gcf, 'resize', 'off')
% set(gcf, 'units', 'inches', 'position', [440   355   600   450]/user_screen_resolution);
% set(gcf, 'resize', 'off')
% drawnow();

% lgd.Position(1:2) = [0.125, 0.925];
% drawnow();
% subplots = extract_subplots(gcf);
% subplots(2).Position(2) = 0.55;
% subplots(4).Position(2) = 0.55;

% drawnow();

% if save_figures
%     save_figure_to_png_svg_fig(figure_save_directory, 'flow-props', [1, 0, 0]);
% end

% return;

% % ==========================
% % plot histogram of valid trials from each dataset
% % ==========================
% dataset_count = reshape(num_valid_trials_datasets, num_datasets, num_window_resolution)';
% X = categorical(dataset_name_array_plot);
% X = reordercats(X, dataset_name_array_plot);

% figure
% bar(X, dataset_count', 'stacked');
% box off
% ylabel('Count')
% legend('WR=64', 'WR=32', 'location', 'northoutside', 'Orientation', 'horizontal')

% if save_figures
%     save_figure_to_png_svg_fig(figure_save_directory, 'datasets-valid-trials', [true, false, false]);
% end

% % ==========================
% % plot pdf of weights
% % ==========================
% % plot_pdf_weights(bins_w, pdf_w_consolidated, individual_method_array, resampling_method_names_plot_short, colors, line_symbols, user_screen_resolution)
% % plot_pdf_weights(bins_w, {pdf_w_consolidated{resampling_method_index_plot, :}}, individual_method_array, resampling_method_names_plot_short, colors, line_symbols, user_screen_resolution)
% plot_pdf_weights(bins_w, {pdf_w_consolidated{resampling_method_index_plot, :}}, individual_method_array, {'Comb'}, colors, line_symbols, user_screen_resolution)

% if save_figures
%     save_figure_to_png_svg_fig(figure_save_directory, 'pdf-weights-2', [1, 0, 0]);
% end

% ==========================
% plot pdfs of errors and uncertainties
% ==========================

unc_indiv_2 = mat2cell(unc_indiv_all_consolidated, num_valid_trials_consolidated*2, ones(1, num_individual_methods));
% unc_comb_2 = mat2cell(unc_comb_all_consolidated, num_valid_trials_consolidated*2, ones(1, num_resampling_methods));
unc_comb_2 = mat2cell(unc_comb_all_consolidated(:, resampling_method_index_plot), num_valid_trials_consolidated*2, ones(1, 1));
violins = make_violin_plot_error_uncertainty_02(err_all_consolidated, unc_indiv_2, unc_comb_2, ...
                                                err_rms_consolidated, unc_indiv_rms_consolidated, unc_comb_rms_consolidated, ...
                                                bins, individual_method_array, {'Comb'}, ...
                                                colors, user_screen_resolution, max_error_plot);
drawnow();

if save_figures
    save_figure_to_png_svg_fig(figure_save_directory, 'violin-2', [1, 0, 0]);
end

return;

% ===================================
%% comparison of true vs estimated error
% ===================================
error_model = 'gaussian';

% ===================================
%% quantile-quantile plot 
% ===================================
% max_error_threshold = 0.2;
figure
% qq_plot_err_true_est(err_all_consolidated, err_est_indiv_consolidated, err_est_comb_consolidated, ...
%                     individual_method_array, resampling_method_names_plot_short, ...
%                     max_error_threshold, 10, colors, symbols, user_screen_resolution);
qq_plot_err_true_est(err_all_consolidated, err_est_indiv_consolidated, err_est_comb_consolidated(:, resampling_method_index_plot), ...
                    individual_method_array, {'Comb'}, ...
                    max_error_plot, 100, colors, symbols, user_screen_resolution);

drawnow();

% save figures
if save_figures
    save_figure_to_png_svg_fig(figure_save_directory, ['qq-error-true-vs-est-' error_model '-3'], [1, 0, 0]);
    % save_figure_to_png_svg_fig(figure_save_directory, ['qq-error-true-vs-est-' error_pdf], [1, 0, 0]);
    % export_fig(fullfile(figure_write_directory, 'qq-error-true-vs-est.png'), '-r600');
end

return;
% % ===================================
% % total variation distance 
% % ===================================
% figure
% % create array of distances
% % Y = [d_err_est_indiv_consolidated, d_err_est_comb_consolidated];
% Y = [d_err_est_indiv_consolidated, d_err_est_comb_consolidated(resampling_method_index_plot)];
% make_histogram_plot(histogram_distance_categories, Y*100, colors);
% drawnow();
% % save figure
% if save_figures
%     save_figure_to_png_svg_fig(figure_save_directory, ['histogram-distance-' error_model], [1, 0, 0]);
%     % save_figure_to_png_svg_fig(figure_save_directory, ['histogram-distance-' error_pdf '-abs'], [1, 0, 0]);
% else
%     pause(0.1);
% end    



% This script calculates the stereo uncertainty of individual and combined
% schemes, by propagation through the stereo measurement chain.

%%
clear
close all
clc

restoredefaultpath;
addpath(genpath('/scratch/shannon/c/aether/Projects/BOS/general-codes/matlab-codes/'));
addpath('../prana-uncertainty-average-dc-new-im-2-save-deform-cs/');
addpath('../general-codes/');
addpath('../stereo_uncertainty_codes_packaged/');

setup_default_settings;
% dbstop if error

% ============================
%% experiment settings
% ============================
% seconds per frame
spf = 0.001;
% No. of frames
num_frames = 50;
% Starting frame
fstart = 24;
% camera numbers
camera_numbers = [1, 3];
num_cameras = numel(camera_numbers);

% ====================================
%% read/write settings
% ====================================
% top level read directory
top_read_directory = '/scratch/shannon/a/aether/Projects/PIV_PTV_Uncertainty_Quantification/Stereo_uncertainty_work/';
% top level write directory
% top_write_directory = '/scratch/shannon/a/aether/Projects/PIV_PTV_Uncertainty_Quantification/uncertainty-combination/analysis/results/stereo-dataset/cam13/With_selfcal/';
top_write_directory = fullfile('/scratch/shannon/a/aether/Projects/PIV_PTV_Uncertainty_Quantification/uncertainty-combination/', '/analysis/results/stereo-dataset/', ['cam' num2str(camera_numbers(1)) num2str(camera_numbers(2))], '/With_selfcal/');

mkdir_c(top_write_directory);
% Load 2d job file
% job_settings = load(fullfile(top_write_directory, 'cam13_VR_prana_fulljob_withselfcal.mat'));
job_settings = load(fullfile(top_write_directory, ['cam' num2str(camera_numbers(1)) num2str(camera_numbers(2)) '_VR_prana_fulljob_withselfcal.mat']));

% extract calibration job file
caljobfile = job_settings.caljobfile;
% extract 2d job file
planarjob = job_settings.planarjob;

% ====================================
%% processing settings
% ====================================
% camera numbers
camera_numbers = [1, 3];
num_cameras = numel(camera_numbers);
% stereo reconstruction type
rectype = 'Willert';
% pass number
pass_number = 3;
% array of uncertainty methods
individual_method_array = {'IM'; 'MC'; 'CS'};
num_individual_methods = numel(individual_method_array);
% name of the uncertainty methods in prana
individual_method_prana_names = {'Uim'; 'MC'; 'Ucs'};
% start frame
start_frame = 25;

% ====================================
%% resampling settings
% ====================================
% number of particles to remove
% num_particles_remove_array = [5, 10, 15];
percentage_particles_remove = 0.25;
% number of resampling trials
num_resampling_trials = 1e3;
% num overall trials
num_trials = 1e3;
% uncertainty flags for the various methods (0 = False, 1 = True)
uncertainty_flags.ppruncertainty = 0;
uncertainty_flags.miuncertainty = 0;
uncertainty_flags.mcuncertainty = 1;

% ============================
%% resampling settings
% ============================
% resampling method names
resampling_method_name_array = {'remove-paired'; 'remove-random'; 'add-random'};
resampling_method_names_plot = {'R-P'; 'R-R'; 'A-R'};
resampling_method_plot_index = 3;
% resampling_method_name_array = {'add-random'};
num_resampling_methods = numel(resampling_method_name_array);
% number of particles to remove
percentage_particles_remove_array = 0.05:0.05:0.25; %0.025:0.025:0.25; 
num_ppr = numel(percentage_particles_remove_array);
% number of resampling trials
num_resampling_trials = 1e2;
% num overall trials
num_trials = 1e3;

% ========================
%% weight settings
% ========================
% number of bins for weights
num_bins_w = 30;
% bins for weights
bins_w = linspace(0, 1, num_bins_w);
% metric type ('rms', 'iqr', 'std')
metric_name = 'iqr';

% ===============================
%% analysis settings
% ===============================
% minimum allowable error (pix.)
error_min = 1e-4;
% maximum allowable error (pix.)
error_max = [0.5, 0.5, 1.5, 1.5];

% number of bins for histogram
num_bins = 30;
% bins for histograms
bins = cell(1, 4);
for component_index = 1:3
    bins{component_index} = linspace(error_min, error_max(component_index), num_bins);
end
bins{4} = bins{3};

% ===============================
% plot settings
% ===============================
% names for the three components
coordinate_names = {'x'; 'y'; 'z'}; % 'All'};
component_names = {'U'; 'V'; 'W'}; % 'All'};
num_components = 3;
% color scheme for the three components
component_colors = {'r'; 'b'; 'g'; 'm'};
method_colors = lines(3);
% line symbols
line_symbols = {'-'; ':'; '-.'};
% symbols
symbols = {'o', '^', 'v'};
% screen resolution
user_screen_resolution = 113;
% number of data points to skip for the qq plot
num_skip = 8; %4;
% save figures? 
save_figures = 1;

% % ====================================
% %% load planar results
% % ====================================
% fprintf('Loading planar results\n');

% results_planar_all = cell(1, 2);
% jobfile_all = cell(1, 2);
% files_im1 = cell(1, 2);
% files_im2 = cell(1, 2);

% % loop through cameras
% for camera_index = 1:num_cameras
%     fprintf('Camera: %d\n', camera_index);
%     % ====================================
%     %% Load data
%     % ====================================
%     % results directory for the current camera
%     current_results_directory = fullfile(top_write_directory, rectype, ['Camera', num2str(camera_index), filesep], 'vectors');
%     % load results for vectors, errors and uncertainties
%     [results_planar_all{camera_index}, num_snapshots] = load_directory_data(current_results_directory, ['VR*pass' num2str(pass_number, '%d') '*.mat']);
%     % extract number of grid point
%     [num_rows, num_cols] = size(results_planar_all{camera_index}{1}.X);
%     num_grid_points = num_rows * num_cols;    
% end

% % ====================================
% % extract job file properties
% % ====================================
% % extract job file
% jobfile = job_settings.job1;

% % ====================================
% %% load stereo results
% % ====================================
% fprintf('Loading stereo results\n');
% directory containing stereo results
stereo_results_directory = fullfile(top_write_directory, rectype, 'Camera1Camera2_3Cfields',filesep);
% load results for vectors, errors and uncertainties
[results_stereo_all, num_snapshots] = load_directory_data(fullfile(stereo_results_directory, 'vectors'), ['piv*pass_' num2str(pass_number, '%d') '*.mat']);
    
% load errors
fprintf('Loading all errors into memory\n');
filename = fullfile(stereo_results_directory, 'errors.mat');
errors_all = load(filename);

% % load uncertainties
% fprintf('Loading all uncertainties into memory\n');
% filename = fullfile(stereo_results_directory, 'uncertainties.mat');
% uncertainties_all = load(filename);

% % ============================
% % Extract calibration details
% % ============================
% % Magnification in mm/pix
% scaling = job_settings.scaling.wil;
% mx = scaling.xscale;
% my = scaling.yscale;
% % Here the results are for camera1camera3 pair
% cam1 = '1';
% cam2 = '3';
% if strcmp(cam1,'1') && strcmp(cam2,'3')
%     mz = scaling.yscale;
% elseif strcmp(cam1,'2') && strcmp(cam2,'4')
%     mz = scaling.xscale;
% end

% ===============================
%% load resampling results
% ===============================
fprintf('loading resampling results\n');
% directory to load results for this case
current_read_directory = fullfile(stereo_results_directory, 'resampling-new', ...
['trials=' num2str(num_trials, '%d') '_nrs' num2str(num_resampling_trials) '-prewindow']);

% load monte carlo results
load(fullfile(current_read_directory, 'monte_carlo_results.mat'));

% ===============================
%% load weights and combined uncertainty
% ===============================
fprintf('loading weights and combined uncertainty\n');
load(fullfile(current_read_directory, 'combined_uncertainties_3c.mat'));

% ===============================
%% load errors and uncertainties statistics
% ===============================
fprintf('loading errors and uncertainty statistics\n');
% load(fullfile(current_read_directory, 'stereo_errors_uncertainties_statistics.mat'));
% load(fullfile(current_read_directory, 'stereo_errors_uncertainties_statistics_02.mat'));
load(fullfile(current_read_directory, 'stereo_errors_uncertainties_statistics_05.mat'));

figure_save_directory = fullfile(current_read_directory, 'figures-newfilt');
mkdir_c(figure_save_directory);
planar_component_names = {'X'; 'Y'};

% % =============================
% %% plot displacement contours
% % =============================
% cmin = -4; cmax = 4;
% clevels = linspace(cmin, cmax, 25);
% xmin = -32; xmax = 18;
% ymin = -25; ymax = 25;

% frame_index = 50;
% X = errors_all.xgrid;
% Y = errors_all.ygrid;
% U = errors_all.Us(:, :, frame_index);
% V = errors_all.Vs(:, :, frame_index);
% W = errors_all.Ws(:, :, frame_index);

% figure
% subplot(1, 3, 1)
% contourf(X, Y, U, clevels, 'edgecolor', 'none')
% colormap('jet')
% % colorbar
% caxis([cmin, cmax])
% set_axes(gca)
% axis([xmin xmax ymin ymax])
% set(gca, 'xtick', [-25, -10, 5]) 
% set(gca, 'ytick', [-20, 0, 20])
% xlabel('X (mm)')
% ylabel('Y (mm)')
% set(gca, 'TickLength', [0.025, 0.025])
% % set(gca, 'xticklabels', [])
% % set(gca, 'yticklabels', [])
% % axis off
% title('U')

% subplot(1, 3, 2)
% contourf(X, Y, V, clevels, 'edgecolor', 'none')
% colormap('jet')
% % colorbar
% caxis([cmin, cmax])
% set_axes(gca)
% axis([xmin xmax ymin ymax])
% set(gca, 'xtick', [-25, -10, 5]) 
% % set(gca, 'xticklabels', [])
% set(gca, 'yticklabels', [])
% xlabel('X (mm)')
% set(gca, 'TickLength', [0.025, 0.025])
% % axis off
% title('V')

% subplot(1, 3, 3)
% contourf(X, Y, W, clevels, 'edgecolor', 'none')
% colormap('jet')
% cb = colorbar;
% caxis([cmin, cmax])
% set_axes(gca)
% axis([xmin xmax ymin ymax])
% set(gca, 'xtick', [-25, -10, 5]) 
% % set(gca, 'xticklabels', [])
% set(gca, 'yticklabels', [])
% xlabel('X (mm)')
% set(gca, 'TickLength', [0.025, 0.025])
% % axis off
% title('W')

% set(gcf, 'resize', 'off');
% set(gcf, 'units', 'inches', 'position', [183   372   829   238]/user_screen_resolution)
% set(gcf, 'resize', 'off');

% cb.Position = [0.925, 0.15, 0.015, 0.75];
% title(cb, '(pix.)')

% drawnow();

% if save_figures
%     save_figure_to_png_svg_fig(figure_save_directory, 'displacement-contour', [1, 0, 0])
% end

% % % =============================
% % %% qq plot of PLANAR individual and resampled uncertainties
% % % =============================
% % for camera_index = 1:num_cameras
% %     for component_index = 1:2
% %         figure
% %         qq_plot_individual_resampled({unc_planar_individual_all_valid{:, component_index, camera_index}}, {unc_planar_resampled_all_valid{:, component_index, camera_index}}, ...
% %                                     individual_method_array, error_max(component_index), num_skip, method_colors, symbols, user_screen_resolution);
% %         drawnow();
% %         sgtitle(['Camera ' num2str(camera_index) ', ' planar_component_names{component_index}], 'fontsize', 14, 'fontweight', 'bold');        
% %         drawnow();
% %     end
% % end

% % ==========================
% % plot variation of ratio with particle % for just one trial
% % ==========================
% trial_index = 100;
% ppr = percentage_particles_remove_array*100;

% % plot_unc_ratio_with_particle_removal(unc_ratio, unc_ratio_fit, trial_index, ppr, ...
% %                                     individual_method_array, resampling_method_names_plot, ['Uncertainty Ratio ' upper(metric_name)], user_screen_resolution);
% % components = {'x'; 'y'; 'z'};
% % plot_unc_ratio_with_particle_removal_stereo(unc_ratio{trial_index}, unc_ratio_fit{trial_index}, ppr, individual_method_array, resampling_method_names_plot, ...
% %                                             upper(metric_name), components, user_screen_resolution)
% plot_unc_ratio_with_particle_removal_stereo(unc_ratio{trial_index}, unc_ratio_fit{trial_index}, ppr, individual_method_array, ...
%                                             {resampling_method_names_plot{resampling_method_plot_index}}, ...
%                                             upper(metric_name), coordinate_names, line_symbols, symbols, user_screen_resolution)
% if save_figures
%     save_figure_to_png_svg_fig(figure_save_directory, [metric_name '-variation'], [1, 0, 0])
% end

% % ========================
% %% pdf of weights
% % ========================
% combination_method_index = 3;
% figure
% for component_index = 1:num_components
%     coordinate_name = coordinate_names{component_index};    
%     subplot(num_components, 1, component_index)
%     for uncertainty_method_index = 1:num_individual_methods                
%         plot(bins_w(1:end-1), pdf_weights.(coordinate_name){combination_method_index, uncertainty_method_index}, ...
%                 line_symbols{uncertainty_method_index}, 'color', [method_colors(1, :)])
%         hold on
%         % axis limits
%         xlim([0 bins_w(end)]);
%         y_lim{combination_method_index, component_index} = ylim;
%         % annotate figures
%         box off
%         title(component_names{component_index});
%         ylabel('PDF')
%         if component_index == num_components
%             xlabel('Weights') 
%         else
%             set(gca, 'xticklabels', []);            
%         end
%     end
%     % pos = get(gca, 'position');
%     % pos(2) = 0.15;
%     % pos(4) = 0.65;
%     % set(gca, 'position', pos)                    
% end
% set_common_limits_subplot(gcf, 'y', 'auto');
% subplot(num_components, 1, num_components)
% lgd = legend('IM', 'MC', 'CS', 'location', 'northoutside', 'orientation', 'horizontal');

% set(gcf, 'resize', 'off');
% % set(gcf, 'units', 'inches', 'position', [254   501   850   350]/user_screen_resolution);
% set(gcf, 'units', 'inches', 'position', [254   300   350   450]/user_screen_resolution)
% set(gcf, 'resize', 'off');
% lgd.Position(1:2) = [0.15, 0.95];

% drawnow();
% % save figures
% if save_figures
%     save_figure_to_png_svg_fig(figure_save_directory, 'pdf-weights', [1, 0, 0]);
%     % export_fig(fullfile(figure_save_directory, 'pdf-weights.png'), '-r600');
% end

% =============================
%% violin plot of error and uncertainty histograms
% =============================
resampling_method_names_plot{resampling_method_plot_index} = 'Comb';
make_violin_plot_error_uncertainty_all_comp_03(err_all_valid, unc_stereo_individual_all_valid, unc_stereo_combined_all_valid, ...
                                                err_rms, unc_indiv_rms', unc_comb_rms', bins, ...
                                                individual_method_array, resampling_method_names_plot, [resampling_method_plot_index], ...
                                                coordinate_names, component_names, method_colors, user_screen_resolution, error_max);

drawnow();

% save figure
if save_figures
    save_figure_to_png_svg_fig(figure_save_directory, 'error-uncertainty-histograms-violin', [1, 0, 0]);
    % export_fig(fullfile(figure_save_directory, ['error-uncertainty-histograms-violin.png']), '-r600'); 
else
    pause(0.1);
end 

return;
% ===================================
%% quantile-quantile plot of true vs estimated error
% ===================================    
resampling_method_names_plot{resampling_method_plot_index} = 'Comb';
% figure
% qq_plot_err_true_est_separate_all_comp(err_all_valid, err_est_valid_individual, err_est_valid_combined, individual_method_array, ...
%                 {resampling_method_names_plot{:}}, component_names, error_max, num_skip, method_colors, symbols, user_screen_resolution)
figure
% qq_plot_err_true_est_separate_all_comp_02(err_all_valid, err_est_valid_individual, {err_est_valid_combined{resampling_method_plot_index, :}}, individual_method_array, ...
%                 {resampling_method_names_plot{resampling_method_plot_index}}, component_names, error_max, num_skip, method_colors, symbols, user_screen_resolution);
qq_plot_err_true_est_separate_all_comp_03(err_all_valid, err_est_valid_individual, err_est_valid_combined, individual_method_array, ...
                resampling_method_names_plot, [resampling_method_plot_index], coordinate_names, component_names, ...
                error_max, num_skip, method_colors, symbols, user_screen_resolution);
% qq_plot_err_true_est_separate_all_comp_03_vertical(err_all_valid, err_est_valid_individual, err_est_valid_combined, individual_method_array, ...
%                 resampling_method_names_plot, [resampling_method_plot_index], coordinate_names, component_names, ...
%                 error_max, num_skip, method_colors, symbols, user_screen_resolution);

drawnow();

% save figures
if save_figures
    save_figure_to_png_svg_fig(figure_save_directory, 'qq-error-true-vs-est-combined-2', [1, 0, 0]);
    % export_fig(fullfile(figure_save_directory, 'qq-error-true-vs-est-combined.png'), '-r600');
end

return;
% % ===================================
% %% quantile-quantile plot of individual and resampled uncertainties
% % ===================================
% figure
% qq_plot_individual_resampled_all_comp(unc_stereo_individual_all_valid, unc_stereo_resampled_all_valid, ...
%                             individual_method_array, error_max, component_names, ... 
%                             num_skip, method_colors, symbols, user_screen_resolution);

% drawnow();
% % save figures
% if save_figures
%     save_figure_to_png_svg_fig(figure_save_directory, 'qq-unc-individual-vs-resampled', [1, 1, 0]);
%     % export_fig(fullfile(figure_save_directory, ['qq-unc-individual-vs-resampled-' component_names{component_index} '.png']), '-r600');
% end


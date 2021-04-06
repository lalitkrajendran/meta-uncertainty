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

% ====================================
%% read/write settings
% ====================================
% top level read directory
top_read_directory = '/scratch/shannon/a/aether/Projects/PIV_PTV_Uncertainty_Quantification/Stereo_uncertainty_work/';
% top level write directory
top_write_directory = '/scratch/shannon/a/aether/Projects/PIV_PTV_Uncertainty_Quantification/uncertainty-combination/analysis/results/stereo-dataset/cam13/With_selfcal/';
mkdir_c(top_write_directory);
% Load 2d job file
job_settings = load(fullfile(top_write_directory, 'cam13_VR_prana_fulljob_withselfcal.mat'));
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
pass_number = 4;
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
num_trials = 1e4;
% uncertainty flags for the various methods (0 = False, 1 = True)
uncertainty_flags.ppruncertainty = 0;
uncertainty_flags.miuncertainty = 0;
uncertainty_flags.mcuncertainty = 1;

% ====================================
% combination method settings
% ====================================
% array of uncertainty combination model to be used
% 'unweighted', 'global-weight-var', 'global-weight-std',
% 'local-weight-var', 'local-weight-std'
% combination_method_array = {'unwt'; 'var'; 'pd-var'; 'entropy'}; %'prob'};
combination_method_array = {'unwt'; 'pd-var'; 'entropy'}; %'prob'};
num_combination_methods = numel(combination_method_array);
combination_method_names = {'Unwt'; 'PD-Var'; 'Entropy'};

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
% bins for histogram of weights
bins_weights = linspace(0, 1, 25);

% ===============================
% plot settings
% ===============================
% names for the three components
component_names = {'U'; 'V'; 'W'; 'All'};
num_components = 3;
% color scheme for the three components
component_colors = {'r'; 'b'; 'g'; 'm'};
method_colors = lines(2);
% line symbols
line_symbols = {'-'; ':'; '-.'};
% symbols
symbols = {'o', '^', 'v'};
% save figures? 
save_figures = 1;
% screen resolution
user_screen_resolution = 113;
% number of data points to skip for the qq plot
num_skip = 2;

% ====================================
%% load planar results
% ====================================
fprintf('Loading planar results\n');

results_planar_all = cell(1, 2);
jobfile_all = cell(1, 2);
files_im1 = cell(1, 2);
files_im2 = cell(1, 2);

% loop through cameras
for camera_index = 1:num_cameras
    fprintf('Camera: %d\n', camera_index);
    % ====================================
    %% Load data
    % ====================================
    % results directory for the current camera
    current_results_directory = fullfile(top_write_directory, rectype, ['Camera', num2str(camera_index), filesep], 'vectors');
    % load results for vectors, errors and uncertainties
    [results_planar_all{camera_index}, num_snapshots] = load_directory_data(current_results_directory, ['VR*pass' num2str(pass_number, '%d') '*.mat']);
    % extract number of grid point
    [num_rows, num_cols] = size(results_planar_all{camera_index}(1).X);
    num_grid_points = num_rows * num_cols;    
end

% ====================================
% extract job file properties
% ====================================
% extract job file
jobfile = job_settings.job1;

% ====================================
%% load stereo results
% ====================================
fprintf('Loading stereo results\n');
% directory containing stereo results
stereo_results_directory = fullfile(top_write_directory, rectype, 'Camera1Camera2_3Cfields',filesep);
% load results for vectors, errors and uncertainties
[results_stereo_all, num_snapshots] = load_directory_data(fullfile(stereo_results_directory, 'vectors'), ['piv*pass_' num2str(pass_number, '%d') '*.mat']);
    
% load errors
fprintf('Loading all errors into memory\n');
filename = fullfile(stereo_results_directory, 'errors.mat');
errors_all = load(filename);

% load uncertainties
fprintf('Loading all uncertainties into memory\n');
filename = fullfile(stereo_results_directory, 'uncertainties.mat');
uncertainties_all = load(filename);

% ============================
% Extract calibration details
% ============================
% Magnification in mm/pix
scaling = job_settings.scaling.wil;
mx = scaling.xscale;
my = scaling.yscale;
% Here the results are for camera1camera3 pair
cam1 = '1';
cam2 = '3';
if strcmp(cam1,'1') && strcmp(cam2,'3')
    mz = scaling.yscale;
elseif strcmp(cam1,'2') && strcmp(cam2,'4')
    mz = scaling.xscale;
end

% ===============================
%% load resampling results
% ===============================
fprintf('loading resampling results\n');
% directory to load results for this case
current_read_directory = fullfile(stereo_results_directory, 'resampling', ...
['trials=' num2str(num_trials, '%d')], ['ppr' num2str(percentage_particles_remove, '%.2f') '_nrs' num2str(num_resampling_trials) '-prewindow']);

% load monte carlo results
load(fullfile(current_read_directory, 'monte_carlo_results.mat'));

% ===============================
%% load weights and combined uncertainty
% ===============================
fprintf('loading weights and combined uncertainty\n');
load(fullfile(current_read_directory, 'combined_uncertainties_3c.mat'));

% % ===============================
% %% load stereo errors and uncertainties
% % ===============================
% fprintf('loading stereo errors and uncertainties\n');
% load(fullfile(current_read_directory, 'stereo_errors_uncertainties.mat'));

% ===============================
%% load errors and uncertainties statistics
% ===============================
fprintf('loading errors and uncertainty statistics\n');
% load(fullfile(current_read_directory, 'stereo_errors_uncertainties_statistics.mat'));
% load(fullfile(current_read_directory, 'stereo_errors_uncertainties_statistics_02.mat'));
load(fullfile(current_read_directory, 'stereo_errors_uncertainties_statistics_03.mat'));

figure_save_directory = fullfile(current_read_directory, 'figures-3cweights-abs');
mkdir_c(figure_save_directory);


% =============================
%% violin plot of error and uncertainty histograms
% =============================
make_violin_plot_error_uncertainty_all_comp(err_all_valid, unc_stereo_individual_all_valid, unc_stereo_combined_all_valid, ...
                                    rms_err, rms_unc_individual, rms_unc_combined, ...
                                    individual_method_array, {combination_method_names{:}}, ...
                                    component_names, method_colors, user_screen_resolution, error_max);

% ===================================
%% quantile-quantile plot of true vs estimated error
% ===================================    
figure
qq_plot_err_true_est_separate_all_comp(err_all_valid, err_est_valid_individual, err_est_valid_combined, individual_method_array, ...
                {combination_method_names{:}}, component_names, error_max, num_skip, method_colors, symbols, user_screen_resolution)
% save figures
if save_figures
    save_figure_to_png_svg_fig(figure_write_directory, 'qq-error-true-vs-est-combined', [1, 1, 1]);
    % export_fig(fullfile(figure_save_directory, 'qq-error-true-vs-est-combined.png'), '-r600');
end

return;
for component_index = 1:(num_components + 1)
    % % =============================
    % %% violin plot of error and uncertainty histograms
    % % =============================
    make_violin_plot_error_uncertainty(abs(err_all_valid{component_index}), {unc_stereo_individual_all_valid{:, component_index}}, {unc_stereo_combined_all_valid{:, component_index}}, ...
                                        rms_err(component_index), rms_unc_individual(:, component_index), rms_unc_combined(:, component_index), ...
                                        individual_method_array, {combination_method_names{:}}, ...
                                        method_colors, user_screen_resolution, error_max(component_index));
    
    % title(component_names{component_index});    
    % drawnow();

    % % save figure
    % if save_figures
    %     % save_figure_to_png_eps_fig(figure_write_directory, 'error-uncertainty-histograms-violin', [1, 0, 0]);
    %     export_fig(fullfile(figure_save_directory, ['error-uncertainty-histograms-violin-' component_names{component_index} '.png']), '-r600'); 
    % else
    %     pause(0.1);
    % end

    % ===================================
    %% quantile-quantile plot of true vs estimated error
    % ===================================    
    figure
    qq_plot_err_true_est_separate(abs(err_all_valid{component_index}), {err_est_valid_individual{:, component_index}}, {err_est_valid_combined{:, component_index}}, ...
                    individual_method_array, {combination_method_names{:}}, error_max(component_index), num_skip, method_colors, symbols, user_screen_resolution);

    % sgtitle(component_names{component_index}, 'fontweight', 'bold');    
    subplot(1, 6, 1)
    text(-error_max(component_index)*0.75, error_max(component_index)/2, component_names{component_index}, 'fontsize', 14, 'fontweight', 'bold')
    drawnow();

    % save figures
    if save_figures
        % save_figure_to_png_eps_fig(figure_write_directory, 'qq-error-true-vs-est', [1, 0, 0]);
        export_fig(fullfile(figure_save_directory, ['qq-error-true-vs-est-' component_names{component_index} '.png']), '-r600');
    end

    % % ===================================
    % %% quantile-quantile plot of individual and resampled uncertainties
    % % ===================================
    % figure
    % qq_plot_individual_resampled({unc_stereo_individual_all_valid{:, component_index}}, {unc_stereo_resampled_all_valid{:, component_index}}, ...
    %                             individual_method_array, error_max(component_index), num_skip, method_colors, symbols, user_screen_resolution);

    % subplot(1, 3, 1)
    % text(-error_max(component_index)/2, error_max(component_index)/2, component_names{component_index}, 'fontsize', 14, 'fontweight', 'bold')
    % % sgtitle(component_names{component_index}, 'fontweight', 'bold');    
    % drawnow();

    % % save figures
    % if save_figures
    %     % save_figure_to_png_eps_fig(figure_write_directory, 'qq-error-true-vs-est', [1, 0, 0]);
    %     export_fig(fullfile(figure_save_directory, ['qq-unc-individual-vs-resampled-' component_names{component_index} '.png']), '-r600');
    % end

end 

% ========================
%% pdf of weights
% ========================
figure
for component_index = 1:num_components
    for combination_method_index = 1:num_combination_methods
        for uncertainty_method_index = 1:num_individual_methods
            subplot_index = (combination_method_index - 1) * 3 + component_index;
            subplot(num_combination_methods, 3, subplot_index)
            if combination_method_index == 1
                plot([1/3, 1/3], [0, 10], line_symbols{uncertainty_method_index}, 'color', component_colors{component_index})
            else
                plot(bins_weights(1:end-1), pdf_weights{combination_method_index, uncertainty_method_index}(component_index, :), ...
                    line_symbols{uncertainty_method_index}, 'color', component_colors{component_index})
            end
            hold on
            % axis limits
            xlim([0 bins_weights(end)]);
            y_lim{combination_method_index, component_index} = ylim;
            % annotate figures
            box off
            if combination_method_index == 1
                title(component_names{component_index});
            end

            if combination_method_index == num_combination_methods
                xlabel('Weights')
            else
                set(gca, 'xticklabel', []);
            end
            
            if component_index == 1 && uncertainty_method_index == num_individual_methods
                % ylabel('PDF');
                text(-0.5, y_lim{combination_method_index, component_index}(2)/2, combination_method_names{combination_method_index}, 'fontsize', 18, 'fontweight', 'bold')
            end    
        end
    end
end

subplot(num_combination_methods, num_components, 1)
lgd = legend('IM', 'MC', 'CS', 'location', 'northoutside', 'orientation', 'horizontal');
lgd.Position(1) = 0.33;
lgd.Position(2) = 0.96;

set(gcf, 'resize', 'off');
drawnow();
% set(gcf, 'units', 'inches', 'position', [310    10   824   800]/user_screen_resolution);
set(gcf, 'units', 'inches', 'position', [315  88  1103  795]/user_screen_resolution);
drawnow();
% save figures
if save_figures
    % save_figure_to_png_eps_fig(figure_write_directory, 'pdf-weights', [1, 0, 0]);
    export_fig(fullfile(figure_save_directory, 'pdf-weights.png'), '-r600');
end

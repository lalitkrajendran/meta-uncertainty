clear
close all
clc

%% set mount directory

if ismac
    mount_directory_a = '/Volumes/aether/';
    mount_directory_c = '/Volumes/aether_c/';
else
    mount_directory_a = '/scratch/shannon/a/aether/';
    mount_directory_c = '/scratch/shannon/c/aether';
end

restoredefaultpath;
addpath(genpath(fullfile(mount_directory_c, 'Projects/BOS/general-codes/matlab-codes/')));
setup_default_settings;

% ========================
%% read/write settings
% ========================

% window resolution
window_resolution_array = [32, 64];
num_window_resolution = numel(window_resolution_array);

% top level directory for this project
top_project_directory = fullfile(mount_directory_a, 'Projects/PIV_PTV_Uncertainty_Quantification/uncertainty-combination/');
% directory containing files to be read
top_read_directory = fullfile(top_project_directory, 'analysis/results/', 'experiment-new');

% array of data set names
dataset_name_array = {'PivChal03B'; 'PivChal05B'; 'stagnation_flow'; 'Vortex_Ring'; 'Jetdata'};
num_datasets = numel(dataset_name_array);

% array of uncertainty methods
uncertainty_method_array = {'IM'; 'MC'; 'CS'};
num_uncertainty_methods = numel(uncertainty_method_array);

% array of uncertainty combination model to be used
% 'unweighted', 'global-weight-var', 'global-weight-std',
% 'local-weight-var', 'local-weight-std'
combination_method_array = {'unwt'; 'var-covar'; 'entropy'}; % 'prob'};
num_combination_methods = numel(combination_method_array);
combination_method_names = {'Unwt'; 'Var-Covar'; 'Entropy'};

% directory where results of this analysis are to be saved
top_write_directory = fullfile(top_project_directory, 'analysis/results/', 'monte-carlo', 'individual-datasets');
mkdir_c(top_write_directory);

% total number of variables to compare (including error)
num_total_methods = 1+num_uncertainty_methods+num_combination_methods;

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

% ========================
%% plot settings
% ========================
% user screen resolution
user_screen_resolution = 113;
% save_figure? (true/false)
save_figures = false;
% symbols
symbols = {'o', '^', 'v'};
% colors
colors = lines(3);
% line symbols
line_symbols = {'-'; ':'; '-.'};
% dataset names to use for plots
dataset_name_array_plot = {'03B'; '05B'; 'SF'; 'VR'; 'Jet'};

figure
% ===============================
%% loop through datasets
% ===============================
for dataset_index = 1:num_datasets
    dataset_name = dataset_name_array{dataset_index};
    fprintf('dataset: %s\n', dataset_name);
    results = cell(1, num_window_resolution);
    for window_resolution_index = 1:num_window_resolution
        fprintf('WS: %d\n', window_resolution_index);
        % ===============================
        %% directory settings for this case
        % ===============================
        % directory to load results for this case
        current_read_directory = fullfile(top_write_directory, dataset_name, ['WS' num2str(window_resolution_index)], ...
        ['trials=' num2str(num_trials, '%d')], ['ppr' num2str(percentage_particles_remove, '%.2f') '_nrs' num2str(num_resampling_trials)], ...
        ['max_error=' num2str(max_error_threshold, '%.2f') 'pix']);

        %% load results
        results{window_resolution_index} = load(fullfile(current_read_directory, 'error_uncertainty_statistics.mat'));
    end
    % =============================
    %% violin plot of error and uncertainty histograms
    % =============================
    ax = subplot(num_datasets, 1, dataset_index);
    hold on
    % set subplot position
    set(ax, 'Position', [0.15, 0.85 - 0.2 * (dataset_index - 1),  0.65, 0.125]);
    % set(ax, 'Position', [0.15, 0.85 - 0.15 * (dataset_index + 1 - 1),  0.65, 0.1]);

    % face colors
    color_all = cell(1, num_total_methods); %numel(violins));
    % plot violins and rms values
    % add lines corresponding to rms
    violins = cell(1, num_total_methods);
    for violin_index = 1:num_total_methods
        % x co-ordinates of current violin plot
        % error
        if violin_index == 1
            violins{violin_index} = violinplot_single_asymm([abs(results{1}.err_all_valid)', abs(results{2}.err_all_valid)'], violin_index, 0.3, [0, 0, 0]);
            y1 = results{1}.err_rms;
            y2 = results{2}.err_rms;
        % individual uncertainty methods
        elseif violin_index > 1 && violin_index <= 1 + num_uncertainty_methods
            violins{violin_index} = violinplot_single_asymm([results{1}.sigma_all_valid{violin_index-1}', results{2}.sigma_all_valid{violin_index-1}'], violin_index, 0.3, colors(1, :));
            plot([violin_index, violin_index], [0, max_error_threshold], 'color', [0, 0, 0, 0.5], 'linewidth', 0.5);
            y1 = results{1}.sigma_rms(violin_index - 1);
            y2 = results{2}.sigma_rms(violin_index - 1);
        % combined uncertainty methods
        elseif violin_index > 4
            violins{violin_index} = violinplot_single_asymm([results{1}.unc_combined_all_valid{violin_index-4}', results{2}.unc_combined_all_valid{violin_index-4}'], violin_index, 0.3, colors(2, :));
            plot([violin_index, violin_index], [0, max_error_threshold], 'color', [0, 0, 0, 0.5], 'linewidth', 0.5);
            y1 = results{1}.sigma_rms_comb(violin_index - 4);
            y2 = results{2}.sigma_rms_comb(violin_index - 4);
        end

        %% adjust violin properties
        color_all{violin_index} = violins{violin_index}.FaceColor;
        violins{violin_index}.FaceAlpha = 0.25;
        violins{violin_index}.EdgeColor = [1, 1, 1];
        
        %% extract plotted points
        x = violins{violin_index}.XData;
        y = violins{violin_index}.YData;

        %% plot dividing line
        plot([violin_index, violin_index], [min(y), max(y)], 'color', color_all{violin_index}, 'linewidth', 0.5);

        %% plot rms value
        % plot([min(x), max(x)], y*[1, 1], 'color', [violins{violin_index}.FaceColor, violins{violin_index}.FaceAlpha+0.2], 'linewidth', 3)
        plot([min(x), violin_index], y1*[1, 1], 'color', [violins{violin_index}.FaceColor, violins{violin_index}.FaceAlpha+0.2], 'linewidth', 2)
        plot([violin_index, max(x)], y2*[1, 1], 'color', [violins{violin_index}.FaceColor, violins{violin_index}.FaceAlpha+0.2], 'linewidth', 2)

        % plot rms error
        plot([min(x), violin_index], results{1}.err_rms*[1, 1], 'color', [0, 0, 0, violins{violin_index}.FaceAlpha+0.2], 'linewidth', 2)
        plot([violin_index, max(x)], results{2}.err_rms*[1, 1], 'color', [0, 0, 0, violins{violin_index}.FaceAlpha+0.2], 'linewidth', 2)
    end

    % adjust limits
    % ylim([0 max_error_threshold])
    ylim([0 0.15])

    pause(0.1);
    % turn off x axis line
    % ax = gca;
    ax.XAxis.Axle.Visible = 'off';
    ax.XAxis.TickLength = [0 0];
    if dataset_index < num_datasets
        ax.XAxis.TickLabel = [];
    end
    % turn off y axis line
    ax.YAxis.Axle.Visible = 'off';
    % ax.YAxis.TickLength = [0, 0];
    ax.YAxis.TickLength = [0.005 0.005];

    % add dataset name
    annotation('textbox', [0.85, 0.9 - 0.2 * (dataset_index - 1), 0, 0], 'string', dataset_name_array_plot{dataset_index}, 'fontsize', 14, 'fontweight', 'bold')    
    % annotation('textbox', [0.85, 0.9 - 0.15 * (dataset_index + 1 - 1), 0, 0], 'string', dataset_name_array_plot{dataset_index}, 'fontsize', 12, 'fontweight', 'bold')    
end

% xlabel
xl = set(gca, 'xticklabel', {'Error'; 'IM'; 'MC'; 'CS'; 'Unwt'; 'Var'; 'Entropy'});
%%
% % adjust figure position
% set(gcf, 'resize', 'off');
% pause(0.1);
% set(gcf, 'units', 'inches', 'Position', [352   526   895   392]/user_screen_resolution)
% set(gca, 'units', 'pix', 'fontsize', 11);

% % save figure
% ['WS' num2str(window_resolution_index)],
% if save_figures
%     % save_figure_to_png_eps_fig(figure_write_directory, 'error-uncertainty-histograms-violin', [1, 0, 0]);
%     export_fig(fullfile(figure_write_directory, 'error-uncertainty-histograms-violin.png'), '-r600'); 
% else
%     pause(0.1);
    % end
% annotate y axis
yl = ylabel('Error/Uncertainty (pix.)', 'fontsize', 16);
yl.Position = [0.2, 0.45, -1];
yl.FontWeight = 'bold';
ax.XAxis.FontSize = 14;
ax.XAxis.FontWeight = 'bold';

set(gcf, 'resize', 'off');
pause(0.1);
set(gcf, 'position', [332    54   855   690])

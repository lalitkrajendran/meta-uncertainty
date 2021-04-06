clear
close all
clc

restoredefaultpath;
addpath(genpath('/scratch/shannon/c/aether/Projects/BOS/general-codes/matlab-codes/'));
setup_default_settings;
addpath('../prana-uncertainty-average-dc-new-im-2-save-deform-cs/');

% ========================
%% read/write settings
% ========================

% window resolution
window_resolution_array = [32, 64];
num_window_resolution = numel(window_resolution_array);

% directory containing images
top_image_directory = fullfile('/scratch/shannon/a/aether/Projects/PIV_PTV_Uncertainty_Quantification/Planar_uncertainty_work/Images/');
% directory containing files to be read
top_read_directory = fullfile('/scratch/shannon/a/aether/Projects/PIV_PTV_Uncertainty_Quantification/uncertainty-combination/analysis/results/', 'experiment-new');

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

% methods to calculate histogram distances
histogram_distance_methods = {'total_variation_distance'; 'chi_square_statistics'; 'kolmogorov_smirnov_distance'; ...
    'hellinger_distance'}; % 'kullback_leibler_divergence'};
num_distance_methods = numel(histogram_distance_methods);

% directory where results of this analysis are to be saved
top_write_directory = fullfile('/scratch/shannon/a/aether/Projects/PIV_PTV_Uncertainty_Quantification/uncertainty-combination/analysis/results/', 'monte-carlo', 'individual-datasets');
mkdir_c(top_write_directory);

% ========================
%% statistical analysis settings
% ========================

% number of trials
num_trials = 1e2;
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

% ===============================
%% loop through window resolutions
% ===============================
for window_resolution_index = 1:num_window_resolution
    fprintf('WS: %d\n', window_resolution_index);
    % figure_write_directory = fullfile(top_write_directory, 'consolidated-figures', ...
    % ['trials=' num2str(num_trials, '%d')], ['ppr' num2str(percentage_particles_remove, '%.2f') '_nrs' num2str(num_resampling_trials)], ...
    % ['max_error=' num2str(max_error_threshold, '%.2f') 'pix'], 'figures');
    % mkdir_c(figure_write_directory);
    figure
    % ===============================
    %% loop through datasets
    % ===============================
    for dataset_index = 1:num_datasets
        dataset_name = dataset_name_array{dataset_index};
        fprintf('dataset: %s\n', dataset_name);
        % ===============================
        %% directory settings for this case
        % ===============================
        % directory to load results for this case
        current_read_directory = fullfile(top_write_directory, dataset_name, ['WS' num2str(window_resolution_index)], ...
        ['trials=' num2str(num_trials, '%d')], ['ppr' num2str(percentage_particles_remove, '%.2f') '_nrs' num2str(num_resampling_trials)], ...
        ['max_error=' num2str(max_error_threshold, '%.2f') 'pix']);

        %% load results
        load(fullfile(current_read_directory, 'error_uncertainty_statistics.mat'));

        % =============================
        %% violin plot of error and uncertainty histograms
        % =============================
        ax = subplot(num_datasets, 1, dataset_index);
        % set subplot position
        set(ax, 'Position', [0.15, 0.85 - 0.2 * (dataset_index - 1),  0.65, 0.125]);
        % make violin plot
        violins = violinplot([abs(err_all_valid)', sigma_all_valid{1}', sigma_all_valid{2}', sigma_all_valid{3}', ...
            unc_combined_all_valid{1}', unc_combined_all_valid{2}', unc_combined_all_valid{3}'], ...
            {'Error'; 'IM'; 'MC'; 'CS'; 'Unwt'; 'Var'; 'Entropy'}, ...
            'showdata', false, 'shownotches', false, 'edgecolor', [1, 1, 1], 'ViolinAlpha', 0.5);

        % face colors
        color_all = cell(1, numel(violins));

        % add lines corresponding to rms
        for violin_index = 1:numel(violins)
            % x co-ordinates of current violin plot
            x = violins(violin_index).ViolinPlot.XData;
            % error
            if violin_index == 1
                violins(violin_index).ViolinColor = [0, 0, 0];
                y = err_rms;
                % x = [violins(violin_index).ViolinPlot.XData; violins(numel(violins)).ViolinPlot.XData];
                plot([min(violins(violin_index).ViolinPlot.XData), max(violins(numel(violins)).ViolinPlot.XData)], y*[1, 1], '--', 'color', [0, 0, 0, 0.2])
            % individual uncertainty methods
            elseif violin_index > 1 && violin_index <= 4
                violins(violin_index).ViolinColor = colors(1, :);
                % violins(violin_index).ViolinColor = colors_blue{violin_index-1};
                y = sigma_rms(violin_index - 1);
            % combined uncertainty methods
            elseif violin_index > 4
                violins(violin_index).ViolinColor = colors(2, :);
                % violins(violin_index).ViolinColor = colors_red{violin_index-4};
                y = sigma_rms_comb(violin_index - 4);
            end
            
            %% adjust violin properties
            color_all{violin_index} = violins(violin_index).ViolinColor;
            violins(violin_index).ViolinAlpha = 0.25;
            violins(violin_index).BoxColor = violins(violin_index).ViolinColor;
            violins(violin_index).BoxWidth = 0.005;
            
            %% plot rms value
            plot([min(x), max(x)], y*[1, 1], 'color', [violins(violin_index).ViolinColor, violins(violin_index).ViolinAlpha+0.2], 'linewidth', 3)
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
    end
    % annotate y axis
    yl = ylabel('Error/Uncertainty (pix.)', 'fontsize', 16);
    yl.Position = [0.2, 0.45, -1];
    set(gcf, 'resize', 'off');
    pause(0.1);
    set(gcf, 'position', [332    54   855   690])
    return;
end

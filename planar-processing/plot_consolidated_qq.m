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
addpath('../general-codes')

% ========================
%% read/write settings
% ========================

% window resolution
window_resolution_array = [32, 64];
num_window_resolution = numel(window_resolution_array);

% top level directory for this project
top_project_directory = fullfile(mount_directory_a, 'Projects/PIV_PTV_Uncertainty_Quantification/uncertainty-combination/');
% directory containing files to be read
top_read_directory = fullfile(top_project_directory, 'analysis/results/planar-dataset', 'experiment-new');

% array of data set names
dataset_name_array = {'PivChal03B'; 'PivChal05B'; 'stagnation_flow'; 'Vortex_Ring'; 'Jetdata'};
num_datasets = numel(dataset_name_array);

% array of uncertainty methods
uncertainty_method_array = {'IM'; 'MC'; 'CS'};
num_uncertainty_methods = numel(uncertainty_method_array);

% array of uncertainty combination model to be used
% 'unweighted', 'global-weight-var', 'global-weight-std',
% 'local-weight-var', 'local-weight-std'
% combination_method_array = {'unwt'; 'var-covar'; 'entropy'}; % 'prob'};
combination_method_array = {'unwt'; 'pd-var'; 'entropy'}; % 'prob'};
num_combination_methods = numel(combination_method_array);
% combination_method_names = {'Unwt'; 'Var-Covar'; 'Entropy'};
combination_method_names = {'Unwt'; 'PD-Var'; 'Entropy'};

% directory where results of this analysis are to be saved
top_write_directory = fullfile(top_project_directory, 'analysis/results/planar-dataset', 'monte-carlo', 'individual-datasets');
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
save_figures = 1;
% symbols
symbols = {'o', '^', 'v'};
% colors
colors = lines(3);
% line symbols
line_symbols = {'-'; ':'; '-.'};
% dataset names to use for plots
dataset_name_array_plot = {'03B'; '05B'; 'SF'; 'VR'; 'Jet'};
% number of data points to skip for the qq plot
num_skip = 10;

figure
ax = cell(1, 10);
% [ha, pos] = tight_subplot(2, 5, [0.001, 0.05], 0.1, [0.1, 0.05]);
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
        ['trials=' num2str(num_trials, '%d')], ['ppr' num2str(percentage_particles_remove, '%.2f') '_nrs' num2str(num_resampling_trials) '-prewindow'], ...
        ['max_error=' num2str(max_error_threshold, '%.2f') 'pix']);

        %% load results
        load(fullfile(current_read_directory, 'error_uncertainty_statistics.mat'));

        % =============================
        %% qq plot of error and uncertainty histograms
        % =============================
        subplot_index = (window_resolution_index - 1) * num_datasets + dataset_index;
        subplot(2, num_datasets, subplot_index);
        ax{subplot_index} = gca;
        % ax = ha((window_resolution_index - 1) * num_datasets + dataset_index);
        % axes(ax{subplot_index});
        plot([0 max_error_threshold], [0 max_error_threshold], 'k');
        hold on

        x = abs(err_all_valid)';
        x = x(isfinite(x));
            
        % individual methods
        qq_sigma = [];
        qq_rms_sigma = nans(1, num_uncertainty_methods);
        % num_skip = 10;
        marker_size = 4;
        legend_string = cell(1, num_uncertainty_methods + num_combination_methods);
        % individual methods
        qq_sigma = [];
        for method_index = 1:num_uncertainty_methods
            y = abs(err_est_valid_individual{method_index})';
            l = make_qq_plot(x, y, num_skip, colors(1, :), symbols{method_index}, marker_size);
            
            qq_sigma = [qq_sigma, l(1)];
        end
        % combined methods
        qq_sigma_comb = [];
        qq_rms_sigma_comb = nans(1, num_combination_methods);
        for method_index = 1:num_combination_methods
            y = abs(err_est_valid_combined{method_index})';
            l = make_qq_plot(x, y, num_skip, colors(2, :), symbols{method_index}, marker_size);
            
            qq_sigma_comb = [qq_sigma_comb, l(1)];
        end

        box off
        axis equal

        axis([0 0.1 0 0.1])
        % add dataset name
        if window_resolution_index == 1
            title(dataset_name_array_plot{dataset_index});
            % ax{subplot_index}.Position(2:4) = [0.5, 0.2, 0.35];
            set(gca, 'Position', [ax{subplot_index}.Position(1)-0.03, 0.5, 0.12, 0.35]);
            if dataset_index == 3
                % add legend
                lgd = legend([qq_sigma, qq_sigma_comb], ...
                {' IM', ' MC', ' CS', ' Unwt', ' Var-Covar', ' Entropy'}, 'location', 'northoutside', 'orientation', 'horizontal');
                lgd.Position = [0.25    0.92    0.5490    0.0416];
                lgd.FontSize = 10;
            end
        elseif window_resolution_index == 2
            % ax{subplot_index}.Position(2:4) = [0.1000    0.2    0.35];
            set(gca, 'Position', [ax{subplot_index}.Position(1)-0.03, 0.1, 0.12, 0.35]);
        end
        
        % ax{subplot_index}.Position(1) = ax{subplot_index}.Position(1) - 0.03; 

        if window_resolution_index == 2 && dataset_index == 1
            xlabel('True Error (pix.)');
            ylabel('Estimated Error (pix.)');
        else
            xlabel('');
            ylabel('');
            ax{subplot_index}.XAxis.TickLabels = '';
            ax{subplot_index}.YAxis.TickLabels = '';
            pause(0.1);
        end        
    end
end

pause(0.1);
set(gcf, 'resize', 'off');
set(gcf, 'Position', [185   257   908   433]);

% set(gcf, 'Position', 'units', 'inches', [185   257   908   433]/user_screen_resolution);
% set(gca, 'units', 'pix', 'fontsize', 10);

figure_write_directory = fullfile(top_write_directory, 'consolidated-figures', ...
    ['trials=' num2str(num_trials, '%d')], ['ppr' num2str(percentage_particles_remove, '%.2f') '_nrs' num2str(num_resampling_trials) '-prewindow'], ...
    ['max_error=' num2str(max_error_threshold, '%.2f') 'pix-new']);
mkdir_c(figure_write_directory);
if save_figures
    % export_fig(fullfile(figure_write_directory, 'qq-true-est-err-all-datasets.png'), '-r600'); 
    save_figure_to_png_svg_fig(figure_write_directory, 'qq-true-est-err-all-datasets', [1, 1, 0]); 
end

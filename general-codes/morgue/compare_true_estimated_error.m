clear
close all
clc

restoredefaultpath;
addpath(genpath('/scratch/shannon/c/aether/Projects/BOS/general-codes/matlab-codes/'));
addpath('../prana-uncertainty-average-dc-new-im-2-save-deform-cs/');

% dbstop if error

% ===================================
%% read/write settings
% ===================================

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
top_write_directory = fullfile('/scratch/shannon/a/aether/Projects/PIV_PTV_Uncertainty_Quantification/uncertainty-combination/analysis/results/', 'monte-carlo', 'new-processing-all-datasets');
mkdir_c(top_write_directory);

% ===================================
%% statistical analysis settings
% ===================================

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

% ===================================
%% resampling settings
% ===================================

% number of particles to remove
% num_particles_remove_array = [5, 10, 15];
percentage_particles_remove = 0.25;
% number of resampling trials
num_resampling_trials = 1e3;

% uncertainty flags for the various methods (0 = False, 1 = True)
uncertainty_flags.ppruncertainty = 0;
uncertainty_flags.miuncertainty = 0;
uncertainty_flags.mcuncertainty = 1;

% ===================================
%% plot settings
% ===================================

% save_figure? (true/false)
save_figures = false;
% symbols
symbols = {'o', '^', 'v'};
% colors
colors = lines(3);
% line symbols
line_symbols = {'-'; ':'; '-.'};

%% directory settings for this case

% directory to save results for this case
current_read_directory = fullfile(top_write_directory, ['trials=' num2str(num_trials, '%d')], ['ppr' num2str(percentage_particles_remove, '%.2f') '_nrs' num2str(num_resampling_trials)], 'new', ['max_error=' num2str(max_error_threshold, '%.2f') 'pix']);

figure_write_directory = fullfile(current_read_directory, 'figures');
mkdir_c(figure_write_directory);

%% load results
load(fullfile(current_read_directory, 'error_uncertainty_statistics.mat'));

% ===================================
%% estimate error from uncertainties
% ===================================
% seed random number generator
rng(100);
% individual methods
num_methods = num_uncertainty_methods;

err_est_valid = cell(1, num_methods);
N_err_est = nans(num_methods, numel(bins)-1);
cdf_err_est = nans(num_methods, numel(bins)-1);

for method_index = 1:num_methods
	num_valid_trials = numel(sigma_all_valid{method_index});
	err_est_current = randn(1, num_valid_trials) .* sigma_all_valid{method_index};
	% only retain valid indices
	err_est_valid{method_index} = err_est_current(abs(err_est_current) > min_error_threshold & ...
								abs(err_est_current) < max_error_threshold);	

	% calculate pdf of estimated error
	[N_err_est(method_index, :), ~] = histcounts(abs(err_est_valid{method_index}), bins, 'normalization', 'pdf');
	cdf_err_est(method_index, :) = histcounts(abs(err_est_valid{method_index}), bins, 'normalization', 'cdf');
end


% combined methods
num_methods = num_combination_methods;

err_est_combined_valid = cell(1, num_methods);
N_err_est_combined = nans(num_methods, numel(bins)-1);
cdf_err_est_combined = nans(num_methods, numel(bins)-1);

for method_index = 1:num_methods
	num_valid_trials = numel(unc_combined_all_valid{method_index});
	err_est_current = randn(1, num_valid_trials) .* unc_combined_all_valid{method_index};
	% only retain valid indices
	err_est_combined_valid{method_index} = err_est_current(abs(err_est_current) > min_error_threshold & ...
								abs(err_est_current) < max_error_threshold);	

	% calculate pdf of estimated error
	[N_err_est_combined(method_index, :), ~] = histcounts(abs(err_est_combined_valid{method_index}), bins, 'normalization', 'pdf');
	cdf_err_est_combined(method_index, :) = histcounts(abs(err_est_combined_valid{method_index}), bins, 'normalization', 'cdf');
end

% ===================================
% plot quantile-quantile plot
% ===================================
figure
plot([0 max_error_threshold], [0 max_error_threshold], 'k');
hold on

x = abs(err_all_valid)';
x = x(isfinite(x));
    
% individual methods
qq_sigma = [];
qq_rms_sigma = nans(1, num_uncertainty_methods);
num_skip = 10;
marker_resolution = 4;
legend_string = cell(1, num_uncertainty_methods + num_combination_methods);
for method_index = 1:num_uncertainty_methods
	y = abs(err_est_valid{method_index})';
 
    % make plot
    l = qqplot(x(1:num_skip:end), y(1:num_skip:end));
    % adjust plot
    % set(l(1), 'marker', 'none'); 
    % set(l(1), 'linestyle', line_symbols{method_index});
    % set(l(1), 'color', colors(1,:));

    set(l(1), 'marker', symbols{method_index});
    set(l(1), 'markeredgecolor', colors(1, :));
    set(l(1), 'markerresolution', marker_resolution);
    
    set(l(2), 'visible', 'off');
    set(l(3), 'visible', 'off');
    
    qq_sigma = [qq_sigma, l(1)];
    % calculate RMS deviation
    qq_rms_sigma(method_index) = rms(l(1).YData - l(1).XData, 'omitnan');
    % construct legend string
    legend_string{method_index} = [uncertainty_method_array{method_index} ' = ' num2str(qq_rms_sigma(method_index), '%.3f')];
end

% combined methods
qq_sigma_comb = [];
qq_rms_sigma_comb = nans(1, num_combination_methods);
for method_index = 1:num_combination_methods
    y = abs(err_est_combined_valid{method_index})';
    
    % make plot
    l = qqplot(x(1:num_skip:end), y(1:num_skip:end));
    % adjust plot
    % set(l(1), 'marker', 'none'); %symbols{method_index});
    % set(l(1), 'linestyle', line_symbols{method_index});
    % set(l(1), 'color', colors(2,:));

    set(l(1), 'marker', symbols{method_index});
    set(l(1), 'markeredgecolor', colors(2, :));
    set(l(1), 'markerresolution', marker_resolution);
    
    set(l(2), 'visible', 'off');
    set(l(3), 'visible', 'off');    
    
    qq_sigma_comb = [qq_sigma_comb, l(1)];
    
    % calculate RMS deviation
    qq_rms_sigma_comb(method_index) = rms(l(1).YData - l(1).XData, 'omitnan');
    % construct legend string
    legend_string{num_uncertainty_methods + method_index} = [convert_string_to_sentence_case(combination_method_array{method_index}) ' = ' num2str(qq_rms_sigma_comb(method_index), '%.3f')];
end

box off
axis equal
axis([0 max_error_threshold 0 max_error_threshold])
xlabel('True Error (pix.)', 'fontresolution', 16)
ylabel('Estimated Error (pix.)', 'fontresolution', 16)
title('Quantile-Quantile Plot', 'fontresolution', 16)
legend([qq_sigma, qq_sigma_comb], ...
    {' IM', ' MC', ' CS', ' Unwt', ' Var-Covar', ' Entropy'}, 'location', 'northoutside', 'NumColumns', 2)
% legend([qq_sigma, qq_sigma_comb], ...
%     legend_string, 'location', 'northoutside', 'NumColumns', 2)

% adjust figure resolution
set(gcf, 'reresolution', 'off');
pause(0.1);
set(gcf, 'Position', [500   423   618   518]);
return;
%%
if save_figures
    save_figure_to_png_eps_fig(figure_write_directory, 'qq-error-true-vs-est', [1, 0, 0]);
end

% ===================================
% plot pdf distributions
% ===================================
figure
ax1 = subplot(2, 1, 1);
for uncertainty_method_index = 1:num_uncertainty_methods    
    plot(bins(1:end-1), N_err_est(uncertainty_method_index, :), line_symbols{uncertainty_method_index}, 'color', colors(1, :))
    hold on
end
plot(bins(1:end-1), N_err, 'k')
xlim([0 max_error_threshold])
% ylim([0 1])
box off
set(gca, 'xticklabel', '');
legend({uncertainty_method_array{:}, 'Error'}, 'location', 'southeast')

ax2 = subplot(2, 1, 2);
for combination_method_index = 1:num_combination_methods  
    plot(bins(1:end-1), N_err_est_combined(combination_method_index, :), line_symbols{combination_method_index}, 'color', colors(2, :))
    hold on
end
plot(bins(1:end-1), N_err, 'k')
xlim([0 max_error_threshold])
% ylim([0 1])
box off
xlabel('(pix.)', 'fontresolution', 16)
legend({combination_method_names{:}, 'Error'}, 'location', 'southeast')

annotation('textarrow',[0.5 0.5],[0.5 0.5],'string','Probability Density Function (PDF)', ...
'HeadStyle','none','LineStyle', 'none', 'TextRotation',90,'Position',[.03 .85 0 0],'Fontresolution',16); %,'FontWeight','bold');
set(gcf, 'reresolution', 'off');
set(gcf, 'position', [680   402   500   600]);

if save_figures
    save_figure_to_png_eps_fig(figure_write_directory, 'pdf-error-true-vs-est', [1, 0, 0]);
end

% ===================================
% plot cdf distributions
% ===================================
figure
ax1 = subplot(2, 1, 1);
for uncertainty_method_index = 1:num_uncertainty_methods    
    plot(bins(1:end-1), cdf_err_est(uncertainty_method_index, :), line_symbols{uncertainty_method_index}, 'color', colors(1, :))
    hold on
end
plot(bins(1:end-1), cdf_err, 'k')
xlim([0 max_error_threshold])
ylim([0 1])
box off
set(gca, 'xticklabel', '');
legend({uncertainty_method_array{:}, 'Error'}, 'location', 'southeast')

ax2 = subplot(2, 1, 2);
for combination_method_index = 1:num_combination_methods  
    plot(bins(1:end-1), cdf_err_est_combined(combination_method_index, :), line_symbols{combination_method_index}, 'color', colors(2, :))
    hold on
end
plot(bins(1:end-1), cdf_err, 'k')
xlim([0 max_error_threshold])
ylim([0 1])
box off
xlabel('(pix.)', 'fontresolution', 16)
legend({combination_method_names{:}, 'Error'}, 'location', 'southeast')

annotation('textarrow',[0.5 0.5],[0.5 0.5],'string','Cumulative Density Function (CDF)', ...
'HeadStyle','none','LineStyle', 'none', 'TextRotation',90,'Position',[.03 .85 0 0],'Fontresolution',16); %,'FontWeight','bold');
set(gcf, 'reresolution', 'off');
set(gcf, 'position', [680   402   500   600]);

if save_figures
    save_figure_to_png_eps_fig(figure_write_directory, 'cdf-error-true-vs-est', [1, 0, 0]);
end

clear
close all
clc

restoredefaultpath;
addpath(genpath('/scratch/shannon/c/aether/Projects/BOS/general-codes/matlab-codes/'));
addpath('../prana-uncertainty-average-dc-new-im-2-save-deform-cs/');

% dbstop if error

%% read/write settings

% window resolution
window_resolution_array = [32, 64];
num_window_size = numel(window_resolution_array);

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

%% statistical analysis settings

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

%% resampling settings

% number of particles to remove
% num_particles_remove_array = [5, 10, 15];
percentage_particles_remove = 0.25;
% number of resampling trials
num_resampling_trials = 1e3;

% uncertainty flags for the various methods (0 = False, 1 = True)
uncertainty_flags.ppruncertainty = 0;
uncertainty_flags.miuncertainty = 0;
uncertainty_flags.mcuncertainty = 1;

%% plot settings

% save_figure? (true/false)
save_figures = true;
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

%% violin plot of error and uncertainty histograms

figure
% make violin plot
% violins = violinplot([abs(err_all_valid)', sigma_all_valid{1}', sigma_all_valid{2}', sigma_all_valid{3}', ...
%     unc_combined_all_valid{1}', unc_combined_all_valid{2}', unc_combined_all_valid{3}', unc_combined_all_valid{4}'], ...
%     {'Error'; 'IM'; 'MC'; 'CS'; 'Unwt'; 'Var'; 'Entropy'; 'Probability'}, ...
%     'showdata', false, 'shownotches', false, 'edgecolor', [1, 1, 1], 'ViolinAlpha', 0.5);

violins = violinplot([abs(err_all_valid)', sigma_all_valid{1}', sigma_all_valid{2}', sigma_all_valid{3}', ...
    unc_combined_all_valid{1}', unc_combined_all_valid{2}', unc_combined_all_valid{3}'], ...
    {'Error'; 'IM'; 'MC'; 'CS'; 'Unwt'; 'Var'; 'Entropy'}, ...
    'showdata', false, 'shownotches', false, 'edgecolor', [1, 1, 1], 'ViolinAlpha', 0.5);

% face colors
color_all = cell(1, numel(violins));
% colors_blue = {[0    0.4470    0.7410]; [0.3000    0.4470    0.7410]; [0.5000    0.4470    0.7410]};
% colors_red = {[1    0.2    0]; [0.9   0.5    0.25]; [1   0.5    0]; [0.800    0.4    0.4]};

% add lines corresponding to rms
for violin_index = 1:numel(violins)
    % x co-ordinates of current violin plot
    x = violins(violin_index).ViolinPlot.XData;
    %% extract rms values
    % error
    if violin_index == 1
        violins(violin_index).ViolinColor = [0, 0, 0];
        y = err_rms;
%         x = [violins(violin_index).ViolinPlot.XData; violins(numel(violins)).ViolinPlot.XData];
        plot([min(violins(violin_index).ViolinPlot.XData), max(violins(numel(violins)).ViolinPlot.XData)], y*[1, 1], '--', 'color', [0, 0, 0, 0.2])
    % individual uncertainty methods
    elseif violin_index > 1 && violin_index <= 4
        violins(violin_index).ViolinColor = colors(1, :);
%         violins(violin_index).ViolinColor = colors_blue{violin_index-1};
        y = sigma_rms(violin_index - 1);
    % combined uncertainty methods
    elseif violin_index > 4
        violins(violin_index).ViolinColor = colors(2, :);
%         violins(violin_index).ViolinColor = colors_red{violin_index-4};
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

% annotate y axis
ylabel('Error/Uncertainty (pix.)', 'fontsize', 16)
% adjust limits
ylim([0 max_error_threshold])

pause(0.1);
% turn off x axis line
ax = gca;
ax.XAxis.Axle.Visible = 'off';
ax.XAxis.TickLength = [0 0];
% turn off y axis line
ax.YAxis.Axle.Visible = 'off';
% ax.YAxis.TickLength = [0, 0];
ax.YAxis.TickLength = [0.005 0.005];

% adjust figure position
set(gcf, 'resize', 'off');
pause(0.1);
set(gcf, 'Position', [352   526   895   392])
% save figure
if save_figures
    save_figure_to_png_eps_fig(figure_write_directory, 'error-uncertainty-histograms-violin', [1, 0, 0]);
else
    pause(0.1);
end
    
%% qqplot of error vs uncertainty

figure
plot([0 max_error_threshold], [0 max_error_threshold], 'k');
hold on

% individual methods
qq_sigma = [];
qq_rms_sigma = nans(1, num_uncertainty_methods);
num_skip = 10;
marker_size = 4;
legend_string = cell(1, num_uncertainty_methods + num_combination_methods);
for method_index = 1:num_uncertainty_methods
    x = abs(err_all_valid)';
    y = sigma_all_valid{method_index}';
    
    % make plot
    l = qqplot(x(1:num_skip:end), y(1:num_skip:end));
    % adjust plot
    set(l(1), 'marker', 'none'); 
    set(l(1), 'linestyle', line_symbols{method_index});
    set(l(1), 'color', colors(1,:));

    % set(l(1), 'marker', symbols{method_index});
    % set(l(1), 'markeredgecolor', colors(1, :));
    % set(l(1), 'markersize', marker_size);
    
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
    x = abs(err_all_valid)';
    y = unc_combined_all_valid{method_index}';
    
    % make plot
    l = qqplot(x(1:num_skip:end), y(1:num_skip:end));
    % adjust plot
    set(l(1), 'marker', 'none'); %symbols{method_index});
    set(l(1), 'linestyle', line_symbols{method_index});
    set(l(1), 'color', colors(2,:));

    % set(l(1), 'marker', symbols{method_index});
    % set(l(1), 'markeredgecolor', colors(2, :));
    % set(l(1), 'markersize', marker_size);
    
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
xlabel('Error (pix.)', 'fontsize', 16)
ylabel('Uncertainty (pix.)', 'fontsize', 16)
title('Quantile-Quantile Plot', 'fontsize', 16)
legend([qq_sigma, qq_sigma_comb], ...
    {' IM', ' MC', ' CS', ' Unwt', ' Var-Covar', ' Entropy'}, 'location', 'northoutside', 'NumColumns', 2)
% legend([qq_sigma, qq_sigma_comb], ...
%     legend_string, 'location', 'northoutside', 'NumColumns', 2)

% adjust figure size
set(gcf, 'resize', 'off');
pause(0.1);
set(gcf, 'Position', [716   423   618   518]);
return;
%%
if save_figures
    save_figure_to_png_eps_fig(figure_write_directory, 'qq-error-uncertainty', [1, 0, 0]);
end

%% display qq plot results

fprintf('==============================\n');
fprintf('RMS from qq plot fits\n');
fprintf('==============================\n');
fprintf('Indiv: %.3f, %.3f, %.3f\n', qq_rms_sigma);    
fprintf('Comb: %.3f, %.3f, %.3f\n', qq_rms_sigma_comb);    

%% pdf of weights

figure
for combination_method_index = 1:num_combination_methods
    ax(combination_method_index) = subplot(num_combination_methods, 1, combination_method_index);
    set(ax(combination_method_index), 'Position', [0.15, 0.7 - 0.3 * (combination_method_index - 1),  0.65, 0.25]);
    % plot histograms for each uncertainty method
    for uncertainty_method_index = 1:num_uncertainty_methods
        plot(bins_weights(1:end-1), N_w{combination_method_index,  uncertainty_method_index}, line_symbols{uncertainty_method_index}, 'color', colors(1, :));
        hold on
    end
    
%     set(gca, 'Position', [0.13, 0.7093, 0.5750, 0.2157])
    box off
    if combination_method_index < num_combination_methods
        set(gca, 'xticklabel', []);
    else
        xl = xlabel('Weights');        
    end
    xlim([0 1])
    
    % get y limits
    y_lim = get(ax(combination_method_index), 'ylim');
    
    plot([1/3, 1/3], y_lim, '--', 'color', [0, 0, 0, 0.5])
%     title(combination_method_names{combination_method_index});    
%     annotation('textbox', [0.6, 0.85 - 0.3 * (combination_method_index - 1), 0, 0], 'string', combination_method_names{combination_method_index}, 'fontsize', 16)
    method_name = convert_string_to_sentence_case(combination_method_array{combination_method_index});
    annotation('textbox', [0.8, 0.85 - 0.3 * (combination_method_index - 1), 0, 0], 'string', method_name, 'fontsize', 14, 'fontweight', 'bold')    
end

axes(ax(combination_method_index))
lgd = legend(uncertainty_method_array, 'location', 'northoutside', 'orientation', 'horizontal', 'fontsize', 14);

lgd.Position = [0.3572 0.9573 0.3199 0.0284];
xl.FontSize = 16;
xl.Position = [0.4742  -1   -1.0000];

annotation('textarrow',[0.5 0.5],[0.5 0.5],'string','Probability Density Function (PDF)', ...
'HeadStyle','none','LineStyle', 'none', 'TextRotation', 90, 'Position', [.04 0.9 0 0], 'FontSize', 16); %, 'FontWeight', 'bold');

set(gcf, 'resize', 'off');
% set(gcf, 'position', [680   370   640   600]);
set(gcf, 'position', [680   370   530   480]);

if save_figures
    save_figure_to_png_eps_fig(figure_write_directory, 'pdf-weights', [1, 0, 0]);
end

%% cdf of error and uncertainty schemes

figure
ax1 = subplot(2, 1, 1);
for uncertainty_method_index = 1:num_uncertainty_methods    
    plot(bins(1:end-1), cdf_sigma(uncertainty_method_index, :), line_symbols{uncertainty_method_index}, 'color', colors(1, :))
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
    plot(bins(1:end-1), cdf_sigma_comb(combination_method_index, :), line_symbols{combination_method_index}, 'color', colors(2, :))
    hold on
end
plot(bins(1:end-1), cdf_err, 'k')
xlim([0 max_error_threshold])
ylim([0 1])
box off
xlabel('(pix.)', 'fontsize', 16)
legend({combination_method_names{:}, 'Error'}, 'location', 'southeast')

annotation('textarrow',[0.5 0.5],[0.5 0.5],'string','Cumulative Density Function (CDF)', ...
'HeadStyle','none','LineStyle', 'none', 'TextRotation',90,'Position',[.03 .85 0 0],'FontSize',16); %,'FontWeight','bold');
set(gcf, 'resize', 'off');
set(gcf, 'position', [680   402   500   600]);

if save_figures
    save_figure_to_png_eps_fig(figure_write_directory, 'cdf-error-uncertainty', [1, 0, 0]);
end

%% pdf of error and uncertainty schemes

figure
ax1 = subplot(2, 1, 1);
for uncertainty_method_index = 1:num_uncertainty_methods    
    plot(bins(1:end-1), N_sigma(uncertainty_method_index, :), line_symbols{uncertainty_method_index}, 'color', colors(1, :))
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
    plot(bins(1:end-1), N_sigma_comb(combination_method_index, :), line_symbols{combination_method_index}, 'color', colors(2, :))
    hold on
end
plot(bins(1:end-1), N_err, 'k')
xlim([0 max_error_threshold])
% ylim([0 1])
box off
xlabel('(pix.)', 'fontsize', 16)
legend({combination_method_names{:}, 'Error'}, 'location', 'southeast')

annotation('textarrow',[0.5 0.5],[0.5 0.5],'string','Probability Density Function (PDF)', ...
'HeadStyle','none','LineStyle', 'none', 'TextRotation',90,'Position',[.03 .85 0 0],'FontSize',16); %,'FontWeight','bold');
set(gcf, 'resize', 'off');
set(gcf, 'position', [680   402   500   600]);

%%
if save_figures
    save_figure_to_png_eps_fig(figure_write_directory, 'pdf-error-uncertainty', [1, 0, 0]);
end

%%
return;
%% coverage plot

X = categorical({'Error', 'IM', 'MC', 'CS', 'Unwt', 'Var-Covar', 'Entropy'}); %, 'Probability'});
% sort in the desired order
X = reordercats(X,{'Error', 'IM', 'MC', 'CS', 'Unwt', 'Var-Covar', 'Entropy'}); %, 'Probability'});
Y = [coverage_err, coverage(:)', coverage_comb(:)'];

figure
% make bar plot
b = bar(X,Y); %,'FaceColor','flat');
hold on
% plot error coverage
l_true = plot(X, coverage_err * ones(1, numel(X)), '--', 'color', [0, 0, 0, 0.2], 'linewidth', 3.0);
% adjust bar properties
b.FaceAlpha = 0.5;
b.BarWidth = 0.5;
b.FaceColor = 'flat';
b.EdgeColor = [1, 1, 1];
b.ShowBaseLine = 'off';
text(1:length(Y),Y,num2str(Y', '%.1f'),'vert','bottom','horiz','center');

% i0 = round(numel(X)/2);
% text(i0,coverage_err+1,num2str(coverage_err, '%.1f'),'vert','bottom','horiz','center');

% adjust color
for i = 1:numel(X)    
    b.CData(i, :) = color_all{i};
end

% set y limit
ylim([0 100])

% turn off box
box off
% annotate y
ylabel('Coverage (%)')
set(gca, 'ycolor', 'none');
% remove ticks
set(gca,'TickLength',[0 0])
ax = gca;
% remove x line and ticks
ax.XAxis.Axle.Visible = 'off';
xtickangle(0)
% figure title
title('Coverage (%)')

% adjust figure size
set(gcf, 'resize', 'off');
pause(0.1);
set(gcf, 'Position', [352   526   895   392])

% save figures
if save_figures
    save_figure_to_png_eps_fig(figure_write_directory, 'coverage-all', [1, 0, 0]);
else
    pause(0.1);
end

%% histogram distance plot

X = categorical({'IM', 'MC', 'CS', 'Unwt', 'Var-Covar', 'Entropy'}); %, 'Probability'});
% sort in the desired order
X = reordercats(X,{'IM', 'MC', 'CS', 'Unwt', 'Var-Covar', 'Entropy'}); %, 'Probability'});

figure
for distance_method_index = 1:num_distance_methods
    % extract distance method name
    current_distance_method = histogram_distance_methods{distance_method_index};
    % create subplot
    subplot(2, 2, distance_method_index)
    % create array of distances
    Y = [d_sigma(:, distance_method_index)', d_sigma_comb(:, distance_method_index)'];
    % find max value
    [maxval, maxloc] = max(Y);
    % find min value
    [minval, minloc] = min(Y);
    
    % normalize distances to the maximum
    Y = Y/maxval * 100;
    % make bar plot
    b = bar(X,Y); %,'FaceColor','flat');
    hold on
    
    % plot error coverage
    plot(X, minval/maxval * 100 * ones(1, numel(X)), '--', 'color', [0, 0, 0, 0.2], 'linewidth', 3.0);

    % adjust bar properties
    b.FaceAlpha = 0.5;
    b.BarWidth = 0.5;
    b.FaceColor = 'flat';
    b.EdgeColor = [1, 1, 1];
    b.ShowBaseLine = 'off';
    % adjust bar color
    for i = 1:numel(X)
        % make max method as black color
        if i == minloc
            b.CData(i, :) = color_all{1};
        else
            b.CData(i, :) = color_all{i+1};
        end
    end
    % add text on top of each bar
    text(1:length(Y),Y,num2str(Y', '%.1f'),'vert','bottom','horiz','center');
	% set axis limit
    ylim([0 100])
    
    % remove box
    box off
    % annotate axis
    ylabel('(%)')
    
    % remove y axis
    set(gca, 'ycolor', 'none');
    set(gca,'TickLength',[0 0])
    
    % turn off x axis
    ax = gca;
    ax.XAxis.Axle.Visible = 'off';
    % remove ticks
    xtickangle(0)
    
    % add title
    t(distance_method_index) = title(convert_string_to_sentence_case(strrep(current_distance_method, '_', ' ')));
    % adjust title position
    t(distance_method_index).Position(2) = 110;


end

% adjust figure size
set(gcf, 'resize', 'off');
pause(0.1);
set(gcf, 'Position', [33         296        1792         725])

% save figure
if save_figures
    save_figure_to_png_eps_fig(figure_write_directory, 'histogram-distance', [1, 0, 0]);
else
    pause(0.1);
end

%% plot rms error vs rms uncertainty

% ==============================
% RMS Error vs RMS Uncertainty
% ==============================
figure
% Error
plot([0 max_error_threshold], [0 max_error_threshold], 'k');
hold on
for method_index = 1:num_uncertainty_methods
    l_sigma(method_index) = plot(err_rms_binwise(method_index, :), sigma_rms_binwise(method_index, :), symbols{method_index}, 'markeredgecolor', colors(1, :));
    plot(err_rms_binwise(method_index, :), sigma_rms_binwise(method_index, :), 'color', [colors(1, :), 0.5]);
end
for method_index = 1:num_combination_methods
    l_sigma_comb(method_index) = plot(err_rms_binwise_comb(method_index, :), sigma_rms_binwise_comb(method_index, :), symbols{method_index}, 'markeredgecolor', colors(2, :));
    plot(err_rms_binwise_comb(method_index, :), sigma_rms_binwise_comb(method_index, :), 'color', [colors(2, :), 0.5]);
end

% grid on
box off
axis equal
axis([0 max_error_threshold 0 max_error_threshold])
xlabel('RMS Error (pix.)')
ylabel('RMS Uncertainty (pix.)')
legend([l_sigma, l_sigma_comb], ...
    {'IM', 'MC', 'CS', 'Unwt', 'Var-Covar', 'Entropy'}, 'location', 'northoutside', 'NumColumns', 2)
% set(gcf, 'Position', [336   511   488   383])

if save_figures
    save_figure_to_png_eps_fig(figure_write_directory, 'rms-error-uncertainty', [true, false, false]);
end


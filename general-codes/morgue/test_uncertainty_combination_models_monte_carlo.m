clear
close all
clc

addpath(genpath('/scratch/shannon/c/aether/Projects/BOS/general-codes/matlab-codes/'));
addpath('prana-uncertainty-average-dc-new-im-2/');
setup_default_settings;
dbstop if error
%% read/write settings

% window resolution
window_resolution_array = [32, 64];
num_window_resolution = numel(window_resolution_array);

% directory containing images
top_image_directory = fullfile('/scratch/shannon/a/aether/Projects/PIV_PTV_Uncertainty_Quantification/Planar_uncertainty_work/Images/');
% directory containing files to be read
top_read_directory = fullfile('/scratch/shannon/a/aether/Projects/PIV_PTV_Uncertainty_Quantification/Planar_uncertainty_work/Results/07_30_2018/');

% array of data set names
dataset_name_array = {'PivChal03B'; 'PivChal05B'; 'stagnation_flow'; 'Vortex_Ring'; 'Jetdata'};
num_datasets = numel(dataset_name_array);

% array of uncertainty combination model to be used
% 'unweighted', 'global-weight-var', 'global-weight-std',
% 'local-weight-var', 'local-weight-std'
uncertainty_combination_method_array = {'unweighted'; 'global-weight-var'; 'local-weight-var'};
num_uncertainty_combination_methods = numel(uncertainty_combination_method_array);

% directory where results of this analysis are to be saved
top_write_directory = fullfile('/scratch/shannon/a/aether/Projects/PIV_PTV_Uncertainty_Quantification/uncertainty-combination/analysis/results/', 'monte-carlo');
if ~exist(top_write_directory, 'dir')
    mkdir(top_write_directory);
end

%% analysis settings

% number of trials
num_trials = 10e3;
% minimum allowable error (pix.)
min_error_threshold = 1e-4;
% maximum allowable error (pix.)
max_error_threshold = 0.2;
% number of bins for histogram
num_bins = 25;
% bins for histograms
bins = linspace(min_error_threshold, max_error_threshold, num_bins);
% number of bins for the coarse histogram
num_bins_coarse = 8;
% save figures? (True/False)
save_figures = true;

%% directory settings for this case

% directory to save results for this case
current_write_directory = fullfile(top_write_directory, ['max_error=' num2str(max_error_threshold, '%.2f') 'pix_trials=' num2str(num_trials, '%d')]);
mkdir_c(current_write_directory);

% directory to save figures
figure_write_directory = fullfile(current_write_directory, 'figures');
mkdir_c(figure_write_directory);

%% pre-load all dataset

fprintf('Loading all datasets into memory\n');

results_all = cell(num_window_resolution, num_datasets);

% loop through window resolutions
for window_resolution_index = 1:num_window_resolution
    window_resolution = window_resolution_array(window_resolution_index);
    % loop through datasets
    for dataset_index = 1:numel(dataset_name_array)
        % name of the data set
        dataset_name = dataset_name_array{dataset_index};
        fprintf('Dataset: %s\n', dataset_name);
        %% Load data

        % load results for vectors, errors and uncertainties
        results_all{window_resolution_index, dataset_index} = load(fullfile(top_read_directory, ['WS' num2str(window_resolution)], '/SOC_gradient/matfiles/', [dataset_name '.mat']));
        % account for peculiarities of each data set
        if dataset_index == 2
            results_all{window_resolution_index, dataset_index}.Xd = fliplr(results_all{window_resolution_index, dataset_index}.Xd');
            results_all{window_resolution_index, dataset_index}.Yd = results_all{window_resolution_index, dataset_index}.Yd';
        end
    end
end


%% allocate memory
fprintf('allocating memory for variables\n');
% dataset address
dataset_index_array = nans(1, num_trials);
window_resolution_index_array = nans(1, num_trials);

% errors
err_Up = nans(1, num_trials);
err_Vp = nans(1, num_trials);
err_Ud = nans(1, num_trials);
err_Vd = nans(1, num_trials);
err_U_comb = nans(1, num_trials);
err_V_comb = nans(1, num_trials);

% uncertainties
sigma_im_x = nans(1, num_trials);
sigma_im_y = nans(1, num_trials);
sigma_mc_x = nans(1, num_trials);
sigma_mc_y = nans(1, num_trials);
sigma_cs_x = nans(1, num_trials);
sigma_cs_y = nans(1, num_trials);

sigma_unwt_x = nans(1, num_trials);
sigma_unwt_y= nans(1, num_trials);
sigma_glob_x = nans(1, num_trials);
sigma_glob_y= nans(1, num_trials);

% weights
w_glob_im_x = nans(1, num_trials);
w_glob_im_y = nans(1, num_trials);
w_glob_mc_x = nans(1, num_trials);
w_glob_mc_y = nans(1, num_trials);
w_glob_cs_x = nans(1, num_trials);
w_glob_cs_y = nans(1, num_trials);

%% loop through trials and perform analysis
fprintf('running monte-carlo\n');

% generate random numbers for overall dataset index
% overall_dataset_index_array = 1 + round(rand(1, num_trials) * (num_datasets * num_window_resolution - 1));
overall_dataset_index_array = 1 + round(rand(1, num_trials) * (num_datasets - 1));

overall_window_resolution_index_array = 1 + round(rand(1, num_trials) * (num_window_resolution - 1));

parfor trial_index = 1:num_trials
    fprintf('trial_index: %d\n', trial_index);

    % calculate random number for data set
    dataset_index_array(trial_index) = randi(num_datasets);
    % calculate random number for window size
    window_resolution_index_array(trial_index) = randi(num_window_resolution);

    % load dataset
    results = results_all{window_resolution_index_array(trial_index), dataset_index_array(trial_index)};
    % calculate size of the results array
    [num_rows, num_cols, num_snapshots] = size(results.Up);

    % calculate random number for the snapshot
    snapshot_index = randi(num_snapshots);

    %% interpolate davis results onto prana grid

    % interpolate davis displacement results onto prana grid
    Ud_interp = interp2(results.Xd, results.Yd, results.Ud(:, :, snapshot_index), results.Xp, results.Yp, 'cubic', 0);
    Vd_interp = interp2(results.Xd, results.Yd, results.Vd(:, :, snapshot_index), results.Xp, results.Yp, 'cubic', 0);

    % interpolate davis displacement uncertainty results onto prana grid
    sigma_Ud_interp = interp2(results.Xd, results.Yd, results.UCSx(:, :, snapshot_index), results.Xp, results.Yp, 'cubic', 0);
    sigma_Vd_interp = interp2(results.Xd, results.Yd, results.UCSy(:, :, snapshot_index), results.Xp, results.Yp, 'cubic', 0);

    % interpolate davis errors on to the prana grid
    err_Ud_interp = interp2(results.Xd, results.Yd, results.err_ud(:, :, snapshot_index), results.Xp, results.Yp, 'cubic', 0);
    err_Vd_interp = interp2(results.Xd, results.Yd, results.err_vd(:, :, snapshot_index), results.Xp, results.Yp, 'cubic', 0);

    %% calculate weights
        
    % assign weights as the inverse of the variance of the
    % uncertainty in field of view
    
    % IM
    w_glob_im_x(trial_index) = 1/var(reshape(results.UIMx(:, :, snapshot_index), num_rows*num_cols, 1), [], 'omitnan');
    w_glob_im_y(trial_index) = 1/var(reshape(results.UIMy(:, :, snapshot_index), num_rows*num_cols, 1), [], 'omitnan');
    % MC
    w_glob_mc_x(trial_index) = 1/var(reshape(results.UMCx(:, :, snapshot_index), num_rows*num_cols, 1), [], 'omitnan');
    w_glob_mc_y(trial_index) = 1/var(reshape(results.UMCy(:, :, snapshot_index), num_rows*num_cols, 1), [], 'omitnan');
    % CS
    w_glob_cs_x(trial_index) = 1/var(sigma_Ud_interp(:), [], 'omitnan');
    w_glob_cs_y(trial_index) = 1/var(sigma_Vd_interp(:), [], 'omitnan');
    %% select random number for grid point in this snapshot
    grid_point_row_index = randi(num_rows);
    grid_point_col_index = randi(num_cols);

    %% extract errors for this grid point
    err_Up(trial_index) = results.err_up(grid_point_row_index, grid_point_col_index, snapshot_index);
    err_Vp(trial_index) = results.err_vp(grid_point_row_index, grid_point_col_index, snapshot_index);

    err_Ud(trial_index) = err_Ud_interp(grid_point_row_index, grid_point_col_index);
    err_Vd(trial_index) = err_Vd_interp(grid_point_row_index, grid_point_col_index);

    %% calculate combined error for this grid point
    err_U_comb(trial_index) = 1/2 * (err_Up(trial_index) + err_Ud(trial_index));
    err_V_comb(trial_index) = 1/2 * (err_Vp(trial_index) + err_Ud(trial_index));

    %% extract individual uncertainties for this grid point
    sigma_im_x(trial_index) = results.UIMx(grid_point_row_index, grid_point_col_index, snapshot_index);
    sigma_im_y(trial_index) = results.UIMy(grid_point_row_index, grid_point_col_index, snapshot_index);

    sigma_mc_x(trial_index) = results.UMCx(grid_point_row_index, grid_point_col_index, snapshot_index);
    sigma_mc_y(trial_index) = results.UMCy(grid_point_row_index, grid_point_col_index, snapshot_index);

    sigma_cs_x(trial_index) = sigma_Ud_interp(grid_point_row_index, grid_point_col_index);
    sigma_cs_y(trial_index) = sigma_Vd_interp(grid_point_row_index, grid_point_col_index);

    %% calculate combined uncertainty for that point for unweighted average
    sigma_unwt_x(trial_index) = 1/3 .* (sigma_im_x(trial_index) + sigma_mc_x(trial_index) + sigma_cs_x(trial_index));
    sigma_unwt_y(trial_index) = 1/3 .* (sigma_im_y(trial_index) + sigma_mc_y(trial_index) + sigma_cs_y(trial_index));
    
    %% calculate combined uncertainty for that point for global weighted average
    sigma_glob_x(trial_index) = 1./(w_glob_im_x(trial_index) + w_glob_mc_x(trial_index) + w_glob_cs_x(trial_index)) .* (w_glob_im_x(trial_index) .* sigma_im_x(trial_index) + w_glob_mc_x(trial_index) .* sigma_mc_x(trial_index) + w_glob_cs_x(trial_index) .* sigma_cs_x(trial_index));
    sigma_glob_y(trial_index) = 1./(w_glob_im_y(trial_index) + w_glob_mc_y(trial_index) + w_glob_cs_y(trial_index)) .* (w_glob_im_y(trial_index) .* sigma_im_y(trial_index) + w_glob_mc_y(trial_index) .* sigma_mc_y(trial_index) + w_glob_cs_y(trial_index) .* sigma_cs_y(trial_index));

end

%% filter out invalid measurements

% remove invalid measurements for the prana error
err_p_all_valid = nan_invalid_measurements([err_Up(:); err_Vp(:)], min_error_threshold, max_error_threshold);
% remove invalid measurements for the davis error
err_d_all_valid = nan_invalid_measurements([err_Ud(:); err_Vd(:)], min_error_threshold, max_error_threshold);
% remove invalid measurements for the combined error
err_comb_all_valid = nan_invalid_measurements([err_U_comb(:); err_V_comb(:)], min_error_threshold, max_error_threshold);

% remove invalid measurements for IM uncertainty
sigma_im_all_valid = nan_invalid_measurements([sigma_im_x(:); sigma_im_y(:)], min_error_threshold, max_error_threshold);
% remove invalid measurements for MC uncertainty
sigma_mc_all_valid = nan_invalid_measurements([sigma_mc_x(:); sigma_mc_y(:)], min_error_threshold, max_error_threshold);
% remove invalid measurements for CS uncertainty
sigma_cs_all_valid = nan_invalid_measurements([sigma_cs_x(:); sigma_cs_y(:)], min_error_threshold, max_error_threshold);
% remove invalid measurements for the unweighted uncertainty
sigma_unwt_all_valid = nan_invalid_measurements([sigma_unwt_x(:); sigma_unwt_y(:)], min_error_threshold, max_error_threshold);
% remove invalid measurements for the globally weighted uncertainty
sigma_glob_all_valid = nan_invalid_measurements([sigma_glob_x(:); sigma_glob_y(:)], min_error_threshold, max_error_threshold);

%% calculate statistics

% calculate rms of prana error
err_p_rms = rms(err_p_all_valid, 'omitnan');
% calculate rms of davis error
err_d_rms = rms(err_d_all_valid, 'omitnan');
% calculate rms of combined error
err_comb_rms = rms(err_comb_all_valid, 'omitnan');

% calculate rms of IM uncertainty
sigma_im_rms = rms(sigma_im_all_valid, 'omitnan');
% calculate rms of MC uncertainty
sigma_mc_rms = rms(sigma_mc_all_valid, 'omitnan');
% calculate rms of CS uncertainty
sigma_cs_rms = rms(sigma_cs_all_valid, 'omitnan');
% calculate rms of unweighted uncertainty
sigma_unwt_rms = rms(sigma_unwt_all_valid, 'omitnan');
% calculate rms of globally weighted uncertainty
sigma_glob_rms = rms(sigma_glob_all_valid, 'omitnan');

%% calculate pdf of the error and uncertainty distributions

% calculate pdf of prana error
[N_err_p, ~] = histcounts(err_p_all_valid, bins, 'normalization', 'pdf');
% calculate pdf of davis error
[N_err_d, ~] = histcounts(err_d_all_valid, bins, 'normalization', 'pdf');
% calculate pdf of combined error
[N_err_comb, ~] = histcounts(err_comb_all_valid, bins, 'normalization', 'pdf');

% calculate pdf of IM uncertainty
[N_sigma_im, ~] = histcounts(sigma_im_all_valid, bins, 'normalization', 'pdf');
% calculate pdf of MC uncertainty
[N_sigma_mc, ~] = histcounts(sigma_mc_all_valid, bins, 'normalization', 'pdf');
% calculate pdf of CS uncertainty
[N_sigma_cs, ~] = histcounts(sigma_cs_all_valid, bins, 'normalization', 'pdf');
% calculate pdf of unweighted uncertainty
[N_sigma_unwt, ~] = histcounts(sigma_unwt_all_valid, bins, 'normalization', 'pdf');
% calculate pdf of globally weighted uncertainty
[N_sigma_glob, ~] = histcounts(sigma_glob_all_valid, bins, 'normalization', 'pdf');

%% calculate coverage

% true prana coverage
cov_p = calculate_coverage_percentage(err_p_all_valid, err_p_rms);
% true davis coverage
cov_d = calculate_coverage_percentage(err_d_all_valid, err_d_rms);
% true combined coverage
cov_comb = calculate_coverage_percentage(err_comb_all_valid, err_comb_rms);

% im coverage
cov_im = calculate_coverage_percentage(err_p_all_valid, sigma_im_all_valid);
% mc coverage
cov_mc = calculate_coverage_percentage(err_p_all_valid, sigma_mc_all_valid);
% cs coverage
cov_cs = calculate_coverage_percentage(err_d_all_valid, sigma_cs_all_valid);
% unweighted coverage
cov_unwt = calculate_coverage_percentage(err_comb_all_valid, sigma_unwt_all_valid);
% globally weighted coverage
cov_glob = calculate_coverage_percentage(err_comb_all_valid, sigma_glob_all_valid);

%% calculate rms error and uncertainty binwise

% calculate binwise rms of error and uncertainty for im
[err_rms_binwise_im, sigma_rms_binwise_im] = calculate_rms_binwise(err_p_all_valid, sigma_im_all_valid, max_error_threshold, num_bins_coarse);
% calculate binwise rms of error and uncertainty for mc
[err_rms_binwise_mc, sigma_rms_binwise_mc] = calculate_rms_binwise(err_p_all_valid, sigma_mc_all_valid, max_error_threshold, num_bins_coarse);
% calculate binwise rms of error and uncertainty for cs
[err_rms_binwise_cs, sigma_rms_binwise_cs] = calculate_rms_binwise(err_d_all_valid, sigma_cs_all_valid, max_error_threshold, num_bins_coarse);
% calculate binwise rms of error and uncertainty for unweighted
[err_rms_binwise_unwt, sigma_rms_binwise_unwt] = calculate_rms_binwise(err_comb_all_valid, sigma_unwt_all_valid, max_error_threshold, num_bins_coarse);
% calculate binwise rms of error and uncertainty for globally weighted
[err_rms_binwise_glob, sigma_rms_binwise_glob] = calculate_rms_binwise(err_comb_all_valid, sigma_glob_all_valid, max_error_threshold, num_bins_coarse);

%% calculate histogram of weights

% bins for histogram
bins_weights = linspace(0, 1, 25);

% calculate normalized weights for all methods, combining x and y weights
w_im_all = [w_glob_im_x./(w_glob_im_x + w_glob_mc_x + w_glob_cs_x), ...
    w_glob_im_y./(w_glob_im_y + w_glob_mc_y + w_glob_cs_y)];
w_mc_all = [w_glob_mc_x./(w_glob_im_x + w_glob_mc_x + w_glob_cs_x), ...
    w_glob_mc_y./(w_glob_im_y + w_glob_mc_y + w_glob_cs_y)];
w_cs_all = [w_glob_cs_x./(w_glob_im_x + w_glob_mc_x + w_glob_cs_x), ...
    w_glob_cs_y./(w_glob_im_y + w_glob_mc_y + w_glob_cs_y)];

% calculate pdf of image matching weights
[N_w_im, ~] = histcounts(w_im_all, bins_weights, 'normalization', 'pdf');
% calculate pdf of mc weights
[N_w_mc, ~] = histcounts(w_mc_all, bins_weights, 'normalization', 'pdf');
% calculate pdf of cs weights
[N_w_cs, ~] = histcounts(w_cs_all, bins_weights, 'normalization', 'pdf');

%% write rms values to file

fprintf('writing statistics to file\n');

fileID = fopen(fullfile(current_write_directory, 'rms-values.txt'), 'w');

fprintf(fileID, '-------------------------------\n');
fprintf(fileID, 'Errors (pix.)\n');
fprintf(fileID, '-------------------------------\n');
fprintf(fileID, 'Prana: %.3f\n', err_p_rms);
fprintf(fileID, 'Davis: %.3f\n', err_d_rms);
fprintf(fileID, 'Combined: %.3f\n', err_comb_rms);

fprintf(fileID, '-------------------------------\n');
fprintf(fileID, 'Uncertainties (pix.)\n');
fprintf(fileID, '-------------------------------\n');

fprintf(fileID, 'IM: %.3f\n', sigma_im_rms);
fprintf(fileID, 'MC: %.3f\n', sigma_mc_rms);
fprintf(fileID, 'CS: %.3f\n', sigma_cs_rms);
fprintf(fileID, 'Unweighted: %.3f\n', sigma_unwt_rms);
fprintf(fileID, 'Global: %.3f\n', sigma_glob_rms);

fclose(fileID);

%% plot histogram of datasets that have been accessed
dataset_count = nans(num_datasets, num_window_resolution);
X = categorical(dataset_name_array);
% find number of trials falling in each dataset
for dataset_index = 1:num_datasets
    for window_resolution_index = 1:num_window_resolution
        indices = find(dataset_index_array == dataset_index);
        dataset_count(dataset_index, window_resolution_index) = sum(window_resolution_index_array(indices) == window_resolution_index);
    end
end

figure
bar(X, dataset_count, 'stacked');
ylabel('Count')
legend('WR=32', 'WR=64', 'location', 'northoutside', 'Orientation', 'horizontal')

if save_figures
    save_figure_to_png_eps_fig(figure_write_directory, 'datasets-accessed', [true, false, false]);
end

%% plot weights
% ====================================
% histogram of normalized weights
% ====================================
figure

bins_plot = bins_weights(1:end-1);
y_max = 10;

% IM
l_im = plot(bins_plot, N_w_im, 'c');
hold on
% MC
l_mc = plot(bins_plot, N_w_mc, 'r');
% CS
l_cs = plot(bins_plot, N_w_cs, 'm');

ylim([0 y_max])
xlabel('Weights')
ylabel('Probability Density')
legend([l_im, l_mc, l_cs], {'IM', 'MC', 'CS'}, 'location', 'northeast')
set(gcf, 'Position', [336   482   583   412])

if save_figures
    save_figure_to_png_eps_fig(figure_write_directory, 'normalized-weights-histogram', [true, false, false]);
end

%% plot histograms of errors and uncertainties with rms 

% =================================
% Error and Uncertainty histograms
% =================================

bins_plot = bins(1:end-1);
y_max = 30;

figure
% average error
l_err_comb = plot(bins_plot, N_err_comb, 'k');
hold on
plot(err_comb_rms*[1 1], [0 y_max], 'k');
% IM
l_im = plot(bins_plot, N_sigma_im, 'c');
plot(sigma_im_rms * [1, 1], [0 y_max], 'c')
% MC
l_mc = plot(bins_plot, N_sigma_mc, 'r');
plot(sigma_mc_rms * [1, 1], [0 y_max], 'r')
% CS
l_cs = plot(bins_plot, N_sigma_cs, 'm');
plot(sigma_cs_rms * [1, 1], [0 y_max], 'm')
% unweighted uncertainty
l_unwt = plot(bins_plot, N_sigma_unwt, 'g');
plot(sigma_unwt_rms*[1 1], [0 y_max], 'g')
% globally weighted uncertainty
l_glob = plot(bins_plot, N_sigma_glob, 'g:');
plot(sigma_glob_rms*[1 1], [0 y_max], 'g:')

ylim([0 y_max])
xlabel('(pix.)')
ylabel('Probability Density')
legend([l_err_comb, l_im, l_mc, l_cs, l_unwt, l_glob], {'e_{comb}', '\sigma_{IM}', '\sigma_{MC}', '\sigma_{CS}', '\sigma_{unwt}', '\sigma_{glob}'}, 'location', 'northeast')
set(gcf, 'Position', [336   482   583   412])

if save_figures
    save_figure_to_png_eps_fig(figure_write_directory, 'histograms', [true, false, false]);
end

%% plot rms error vs rms uncertainty
% ==============================
% RMS Error vs RMS Uncertainty
% ==============================
figure
% Error
plot([0 max_error_threshold], [0 max_error_threshold], 'k--');
hold on
% IM
l_im = plot(err_rms_binwise_im, sigma_rms_binwise_im, 'co-', 'markerfacecolor', 'c');
% MC
l_mc = plot(err_rms_binwise_mc, sigma_rms_binwise_mc, 'ro-', 'markerfacecolor', 'r');
% CS
l_cs = plot(err_rms_binwise_cs, sigma_rms_binwise_cs, 'mo-', 'markerfacecolor', 'm');
% unweighted uncertainty
l_unwt = plot(err_rms_binwise_unwt, sigma_rms_binwise_unwt, 'go-');
% globally weighted uncertainty
l_glob = plot(err_rms_binwise_glob, sigma_rms_binwise_glob, 'g^-');

grid on
axis equal
axis([0 max_error_threshold 0 max_error_threshold])
xlabel('RMS Error (pix.)')
ylabel('RMS Uncertainty (pix.)')
legend([l_im, l_mc, l_cs, l_unwt, l_glob], {'IM', 'MC', 'CS', 'Unwt', 'Glob'}, 'location', 'eastoutside')
set(gcf, 'Position', [336   511   488   383])

if save_figures
    save_figure_to_png_eps_fig(figure_write_directory, 'rms-error-uncertainty', [true, false, false]);
end

%% plot coverage
% ====================================
% coverage
% ====================================

X = categorical({'IM', 'MC', 'CS', 'Unwt', 'Glob'});
% sort in the desired order
X = reordercats(X,{'IM', 'MC', 'CS', 'Unwt', 'Glob'});
Y = [cov_im, cov_mc, cov_cs, cov_unwt, cov_glob];

figure
bar(X,Y); %,'FaceColor','flat');
hold on
l_p = plot(X, cov_p * [1, 1, 1, 1, 1], 'k');
l_d = plot(X, cov_d * [1, 1, 1, 1, 1], 'k--');
l_comb = plot(X, cov_comb * [1, 1, 1, 1, 1], 'k:');

ylim([0 100])
ylabel('Coverage (%)')
legend([l_p, l_d, l_comb], {'Target, Prana'; 'Target, Davis'; 'Target, Comb'}, 'location', 'northoutside')
set(gcf, 'Position', [336   482   583   412])
if save_figures
    save_figure_to_png_eps_fig(figure_write_directory, 'coverage', [true, false, false]);
end

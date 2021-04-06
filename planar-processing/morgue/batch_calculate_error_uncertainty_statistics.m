clear
close all
clc

restoredefaultpath;
addpath ../histogram_distance/
addpath(genpath('/scratch/shannon/c/aether/Projects/BOS/general-codes/matlab-codes/'));
addpath('../prana-uncertainty-average-dc-new-im-2-save-deform-cs/');

% dbstop if error

%% read/write settings

% window resolution
window_resolution_array = [32, 64];
num_window_size = numel(window_resolution_array);

% pass number
pass_number = 4;

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
combination_methods = {'unwt'; 'var'; 'entropy'}; % 'prob'};
num_combination_methods = numel(combination_methods);

% methods to calculate histogram distances
histogram_distance_methods = {'total_variation_distance'; 'chi_square_statistics'; 'kolmogorov_smirnov_distance'; ...
    'hellinger_distance'; 'kullback_leibler_divergence'};
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
bins = linspace(0, max_error_threshold, num_bins);
% number of bins for the coarse histogram
num_bins_coarse = 8;
% bins for histogram of weights
bins_weights = linspace(0, 1, 25);
% small value to replace zero
small_value = 1e-10;
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

%% directory settings for this case

% directory to save results for this case
current_read_directory = fullfile(top_write_directory, ['trials=' num2str(num_trials, '%d')], ['ppr' num2str(percentage_particles_remove, '%.2f') '_nrs' num2str(num_resampling_trials)], 'new');

% directory to save results for this case
current_write_directory = fullfile(current_read_directory, ['max_error=' num2str(max_error_threshold, '%.2f') 'pix']);
mkdir_c(current_write_directory);

%% load results
load(fullfile(current_read_directory, 'monte_carlo_results.mat'));
load(fullfile(current_read_directory, 'combined_uncertainties.mat'));

%% calculate error statistics

% aggregate error
err_all = abs([err_U(:); err_V(:)]);
% remove invalid measurements for the error
err_all_valid = nan_invalid_measurements(err_all, min_error_threshold, max_error_threshold);
err_all_valid = err_all_valid(isfinite(err_all_valid));
% calculate pdf of error
[N_err, ~] = histcounts(err_all_valid , bins, 'normalization', 'pdf');
cdf_err = histcounts(err_all_valid, bins, 'normalization', 'cdf');

% calculate rms of error
err_rms = rms(err_all_valid, 'omitnan');

% coverage of error
coverage_err = calculate_coverage_percentage(err_all_valid, err_rms);

%% calculate uncertainty statistics for individual schemes

sigma_all_valid = cell(1, num_uncertainty_methods);
N_sigma = nans(num_uncertainty_methods, num_bins-1);
cdf_sigma = nans(num_uncertainty_methods, num_bins-1);
sigma_rms = nans(1, num_uncertainty_methods);
coverage = nans(1, num_uncertainty_methods);
err_rms_binwise = nans(num_uncertainty_methods, num_bins_coarse);
sigma_rms_binwise = nans(num_uncertainty_methods, num_bins_coarse);

%% aggregate measurements

sigma_all = {[]; []; []};
for trial_index = 1:num_trials
    if ~isempty(unc_trials{trial_index})
        sigma_all{1} = [sigma_all{1}, unc_trials{trial_index}.imx, unc_trials{trial_index}.imy];
        sigma_all{2} = [sigma_all{2}, unc_trials{trial_index}.mcx, unc_trials{trial_index}.mcy];
        sigma_all{3} = [sigma_all{3}, unc_trials{trial_index}.csx, unc_trials{trial_index}.csy];
    else
        sigma_all{1} = [sigma_all{1}, NaN, NaN];
        sigma_all{2} = [sigma_all{2}, NaN, NaN];
        sigma_all{3} = [sigma_all{3}, NaN, NaN];
    end
end

%% calculate statistics

for uncertainty_method_index = 1:num_uncertainty_methods
    
    % remove invalid measurements
    sigma_all_valid{uncertainty_method_index} = nan_invalid_measurements(sigma_all{uncertainty_method_index}, min_error_threshold, max_error_threshold);
    sigma_all_current = sigma_all_valid{uncertainty_method_index};
    
    % calculate rms of uncertainty
    sigma_rms(uncertainty_method_index) = rms(sigma_all_valid{uncertainty_method_index}, 'omitnan');
    % calculate pdf of IM uncertainty
    N_sigma(uncertainty_method_index, :) = histcounts(sigma_all_current(isfinite(sigma_all_current)), bins, 'normalization', 'pdf');
    cdf_sigma(uncertainty_method_index, :) = histcounts(sigma_all_current(isfinite(sigma_all_current)), bins, 'normalization', 'cdf');
    % calculate coverage
    coverage(uncertainty_method_index) = calculate_coverage_percentage(err_all_valid, sigma_all_valid{uncertainty_method_index});
    % calculate rms error and uncertainty binwise
    [err_rms_binwise(uncertainty_method_index, :), sigma_rms_binwise(uncertainty_method_index, :)] = calculate_rms_binwise(err_all_valid, sigma_all_valid{uncertainty_method_index}, ...
                                            max_error_threshold, num_bins_coarse);

end

%% calculations for combination methods

unc_combined_all_valid = cell(1, num_combination_methods);
N_sigma_comb = nans(num_combination_methods, num_bins-1);
cdf_sigma_comb = nans(num_combination_methods, num_bins-1);
sigma_rms_comb = nans(1, num_combination_methods);
coverage_comb = nans(1, num_combination_methods);
err_rms_binwise_comb = nans(num_combination_methods, num_bins_coarse);
sigma_rms_binwise_comb = nans(num_combination_methods, num_bins_coarse);
N_w = cell(num_combination_methods, num_uncertainty_methods);

%% aggregate measurements
unc_combined_all = {[]; []; []; []};
weights_all = repmat({[]}, 4, 3);
weights_all_valid = weights_all;
for trial_index = 1:num_trials
    unc_combined_current = unc_combined{trial_index};
    weights_current = weights{trial_index};
    
    for combination_method_index = 1:num_combination_methods
        
        % aggregate uncertainties
        if ~isempty(unc_combined_current)
            unc_combined_all{combination_method_index} = [unc_combined_all{combination_method_index}, unc_combined_current{combination_method_index}.x, unc_combined_current{combination_method_index}.y];
        else
            unc_combined_all{combination_method_index} = [unc_combined_all{combination_method_index}, NaN, NaN];
        end
        
        % aggregate weights
        for uncertainty_method_index = 1:num_uncertainty_methods
           if ~isempty(weights_current)
                weights_all{combination_method_index, uncertainty_method_index} = [weights_all{combination_method_index, uncertainty_method_index}, weights_current{combination_method_index}.x(uncertainty_method_index), ...
                                                                    weights_current{combination_method_index}.y(uncertainty_method_index)];
           else
                weights_all{combination_method_index, uncertainty_method_index} = [weights_all{combination_method_index, uncertainty_method_index}, NaN, NaN];
           end
        end
    end
end

%%
for combination_method_index = 1:num_combination_methods
    %% remove invalid measurements for the combined uncertainty
    
    [unc_combined_all_valid{combination_method_index}, valid_indices] = nan_invalid_measurements(unc_combined_all{combination_method_index}, min_error_threshold, max_error_threshold);
    
    %% calculate statistics    

    % calculate rms of unweighted uncertainty
    sigma_rms_comb(combination_method_index) = rms(unc_combined_all_valid{combination_method_index}, 'omitnan');
    
    %% calculate pdf of the error and uncertainty distributions
    sigma_all_current = unc_combined_all_valid{combination_method_index};
    % calculate pdf of unweighted uncertainty
    N_sigma_comb(combination_method_index, :) = histcounts(sigma_all_current(isfinite(sigma_all_current)), bins, 'normalization', 'pdf');
    cdf_sigma_comb(combination_method_index, :) = histcounts(sigma_all_current(isfinite(sigma_all_current)), bins, 'normalization', 'cdf');
    
    %% calculate coverage
    
    % unweighted coverage
    coverage_comb(combination_method_index) = calculate_coverage_percentage(err_all_valid, unc_combined_all_valid{combination_method_index});
    
    %% calculate rms error and uncertainty binwise
    
    % calculate binwise rms of error and uncertainty for unweighted
    [err_rms_binwise_comb(combination_method_index, :), sigma_rms_binwise_comb(combination_method_index, :)] = calculate_rms_binwise(err_all_valid, unc_combined_all_valid{combination_method_index}, max_error_threshold, num_bins_coarse);

    %% calculate histogram of weights
        
    % calculate pdf of weights for the uncertainty schemes
    for uncertainty_method_index = 1:num_uncertainty_methods
        weights_all_valid{combination_method_index, uncertainty_method_index} = weights_all{combination_method_index, uncertainty_method_index}(valid_indices);

        [N_w{combination_method_index, uncertainty_method_index}, ~] = histcounts(weights_all_valid{combination_method_index, uncertainty_method_index}, bins_weights, 'normalization', 'pdf');
%         [N_w{combination_method_index, uncertainty_method_index}, ~] = histcounts(weights_all_valid{combination_method_index, uncertainty_method_index}, bins_weights);
    end
            
end

%% calculate histograms

% scaling factor
scaling_factor = num_trials * (bins(2) - bins(1));

% error
h_err = N_err * scaling_factor;
% replace zero values with small finite number
indices = h_err == 0;
h_err(indices) = small_value;

% individual methods
h_sigma = nans(num_uncertainty_methods, num_bins-1);
% loop through individual methods
for uncertainty_method_index = 1:num_uncertainty_methods
    % calcualte histogram from pdf
    h_sigma(uncertainty_method_index, :) = N_sigma(uncertainty_method_index, :) * scaling_factor;
    % replace zero values with small finite number
    indices = h_sigma(uncertainty_method_index, :) == 0;
    h_sigma(uncertainty_method_index, indices) = small_value;
end

% combined methods
h_sigma_comb = nans(num_combination_methods, num_bins-1);
% loop through combination method
for combination_method_index = 1:num_combination_methods
    % calcualte histogram from pdf
    h_sigma_comb(combination_method_index, :) = N_sigma_comb(combination_method_index, :) * scaling_factor;
    % replace zero values with small finite number
    indices = h_sigma_comb(combination_method_index, :) == 0;
    h_sigma_comb(combination_method_index, indices) = small_value;
end

%% calculate histogram distances

% individual methods
d_sigma = nans(num_uncertainty_methods, num_distance_methods);
% loop through individual methods
for uncertainty_method_index = 1:num_uncertainty_methods
    % loop through distance calculation method
    for distance_method_index = 1:num_distance_methods
        % current method
        current_distance_method = histogram_distance_methods{distance_method_index};
        % calculate distance measure
        d_sigma(uncertainty_method_index, distance_method_index) = pdist2(h_err, h_sigma(uncertainty_method_index, :), str2func(current_distance_method));
    end
end


% combined methods
d_sigma_comb = nans(num_combination_methods, num_distance_methods);
% loop through combination method
for combination_method_index = 1:num_combination_methods
    % loop through distance calculation method
    for distance_method_index = 1:num_distance_methods
        % current method
        current_distance_method = histogram_distance_methods{distance_method_index};
        % calculate distance measure
        d_sigma_comb(combination_method_index, distance_method_index) = pdist2(h_err, h_sigma_comb(combination_method_index, :), str2func(current_distance_method));
    end
end

%% display results

for distance_method_index = 1:num_distance_methods
    fprintf('==============================\n');
    fprintf('%s\n', histogram_distance_methods{distance_method_index});
    fprintf('Indiv: %.2f, %.2f, %.2f\n', d_sigma(:, distance_method_index));    
    fprintf('Comb: %.2f, %.2f, %.2f\n', d_sigma_comb(:, distance_method_index));    
end

%%
filename = fullfile(current_write_directory, 'error_uncertainty_statistics.mat');
save(filename, 'err_all_valid', 'N_err', 'cdf_err', 'err_rms', 'coverage_err', ...
    'sigma_all_valid', 'N_sigma', 'cdf_sigma', 'sigma_rms', 'coverage', 'err_rms_binwise', 'sigma_rms_binwise', ...
    'unc_combined_all_valid', 'weights_all_valid', 'N_sigma_comb', 'cdf_sigma_comb', 'sigma_rms_comb', 'coverage_comb', ...
    'err_rms_binwise_comb', 'sigma_rms_binwise_comb', 'N_w', ...
    'h_err', 'h_sigma', 'h_sigma_comb', 'd_sigma', 'd_sigma_comb');
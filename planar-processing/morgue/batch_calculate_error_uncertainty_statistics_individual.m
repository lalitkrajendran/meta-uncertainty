clear
close all
clc

restoredefaultpath;
addpath('../histogram_distance/') 
addpath(genpath('/scratch/shannon/c/aether/Projects/BOS/general-codes/matlab-codes/'));
addpath('../prana-uncertainty-average-dc-new-im-2-save-deform-cs/');
addpath('../general-codes') 

% ===============================
%% read/write settings
% ===============================
% window resolution
window_resolution_array = [32, 64];
num_window_resolution = numel(window_resolution_array);

% pass number
pass_number = 4;

% directory containing images
top_image_directory = fullfile('/scratch/shannon/a/aether/Projects/PIV_PTV_Uncertainty_Quantification/Planar_uncertainty_work/Images/');
% directory containing files to be read
top_read_directory = fullfile('/scratch/shannon/a/aether/Projects/PIV_PTV_Uncertainty_Quantification/uncertainty-combination/analysis/results/planar-dataset/', 'experiment-new');

% array of data set names
dataset_name_array = {'PivChal03B'; 'PivChal05B'; 'stagnation_flow'; 'Vortex_Ring'; 'Jetdata'};
num_datasets = numel(dataset_name_array);

% array of uncertainty methods
uncertainty_method_array = {'IM'; 'MC'; 'CS'};
num_uncertainty_methods = numel(uncertainty_method_array);

% array of uncertainty combination model to be used
% 'unweighted', 'global-weight-var', 'global-weight-std',
% 'local-weight-var', 'local-weight-std'
combination_method_array = {'unwt'; 'pd-var'; 'entropy'}; % 'prob'};
num_combination_methods = numel(combination_method_array);

% methods to calculate histogram distances
histogram_distance_methods = {'total_variation_distance'; 'chi_square_statistics'; 'kolmogorov_smirnov_distance'; ...
    'hellinger_distance'; 'kullback_leibler_divergence'};
num_distance_methods = numel(histogram_distance_methods);

% directory where results of this analysis are to be saved
top_write_directory = fullfile('/scratch/shannon/a/aether/Projects/PIV_PTV_Uncertainty_Quantification/uncertainty-combination/analysis/results/planar-dataset/', 'monte-carlo', 'individual-datasets');
mkdir_c(top_write_directory);

% ===============================
%% statistical analysis settings
% ===============================
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

% ===============================
%% resampling settings
% ===============================
% number of particles to remove
percentage_particles_remove = 0.25;
% number of resampling trials
num_resampling_trials = 1e3;

% uncertainty flags for the various methods (0 = False, 1 = True)
uncertainty_flags.ppruncertainty = 0;
uncertainty_flags.miuncertainty = 0;
uncertainty_flags.mcuncertainty = 1;

% ===============================
%% declare variables to store aggregated results
% ===============================
err_all_valid_all_datasets = [];
sigma_all_valid_all_datasets = {[]; []; []};
weights_all_valid_all_datasets = {[], [], []; [], [], []; [], [], []};
unc_combined_all_valid_all_datasets = {[]; []; []};
err_est_valid_individual_all_datasets = {[]; []; []};
err_est_valid_combined_all_datasets = {[]; []; []};

% ===============================
%% loop through window resolutions
% ===============================
for window_resolution_index = 1:num_window_resolution
    fprintf('WS: %d\n', window_resolution_index);

    % ===============================
    %% loop through datasets
    % ===============================
    for dataset_index = 1:num_datasets
        fprintf('dataset: %s\n', dataset_name_array{dataset_index});
        % =================================
        %% directory settings for this case
        % =================================
        % directory to load results for this case
        current_read_directory = fullfile(top_write_directory, dataset_name_array{dataset_index}, ['WS' num2str(window_resolution_index)], ...
        ['trials=' num2str(num_trials, '%d')], ['ppr' num2str(percentage_particles_remove, '%.2f') '_nrs' num2str(num_resampling_trials) '-prewindow']);

        % directory to save results for this case
        current_write_directory = fullfile(current_read_directory, ['max_error=' num2str(max_error_threshold, '%.2f') 'pix']);
        mkdir_c(current_write_directory);

        %% load results
        load(fullfile(current_read_directory, 'monte_carlo_results.mat'));
        load(fullfile(current_read_directory, 'combined_uncertainties.mat'));

        % ===============================
        %% calculate error statistics
        % ===============================

        % aggregate error
        err_all = abs([err_U(:); err_V(:)]);
        % remove invalid measurements for the error
        [err_all_valid, valid_indices] = nan_invalid_measurements(err_all, min_error_threshold, max_error_threshold);
        % err_all_valid = err_all_valid(isfinite(err_all_valid));
        % calculate pdf of error
        [pdf_err, ~] = histcounts(err_all_valid(isfinite(err_all_valid)), bins, 'normalization', 'pdf');
        cdf_err = histcounts(err_all_valid(isfinite(err_all_valid)), bins, 'normalization', 'cdf');

        % calculate rms of error
        err_rms = rms(err_all_valid, 'omitnan');

        % coverage of error
        coverage_err = calculate_coverage_percentage(err_all_valid, err_rms);

        % ===============================
        %% calculate uncertainty statistics for individual schemes
        % ===============================
        sigma_all_valid = cell(1, num_uncertainty_methods);
        pdf_sigma = nans(num_uncertainty_methods, num_bins-1);
        cdf_sigma = nans(num_uncertainty_methods, num_bins-1);
        sigma_rms = nans(1, num_uncertainty_methods);
        coverage = nans(1, num_uncertainty_methods);
        err_rms_binwise = nans(num_uncertainty_methods, num_bins_coarse);
        sigma_rms_binwise = nans(num_uncertainty_methods, num_bins_coarse);

        % ===============================
        %% aggregate measurements
        % ===============================
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

        % ===============================
        %% calculate statistics
        % ===============================
        for uncertainty_method_index = 1:num_uncertainty_methods 
            % remove invalid measurements
            sigma_all_valid{uncertainty_method_index} = nan_invalid_measurements(real(sigma_all{uncertainty_method_index}), min_error_threshold, max_error_threshold);
            sigma_all_current = real(sigma_all_valid{uncertainty_method_index});
            
            % calculate rms of uncertainty
            sigma_rms(uncertainty_method_index) = rms(sigma_all_valid{uncertainty_method_index}, 'omitnan');
            % calculate pdf of IM uncertainty
            pdf_sigma(uncertainty_method_index, :) = histcounts(sigma_all_current(isfinite(sigma_all_current)), bins, 'normalization', 'pdf');
            cdf_sigma(uncertainty_method_index, :) = histcounts(sigma_all_current(isfinite(sigma_all_current)), bins, 'normalization', 'cdf');
            % calculate coverage
            coverage(uncertainty_method_index) = calculate_coverage_percentage(err_all_valid, sigma_all_valid{uncertainty_method_index});
            % calculate rms error and uncertainty binwise
            [err_rms_binwise(uncertainty_method_index, :), sigma_rms_binwise(uncertainty_method_index, :)] = calculate_rms_binwise(err_all_valid, sigma_all_valid{uncertainty_method_index}, ...
                                                    max_error_threshold, num_bins_coarse);
        end

        % ===============================
        %% calculations for combination methods
        % ===============================
        unc_combined_all_valid = cell(1, num_combination_methods);
        pdf_sigma_comb = nans(num_combination_methods, num_bins-1);
        cdf_sigma_comb = nans(num_combination_methods, num_bins-1);
        sigma_rms_comb = nans(1, num_combination_methods);
        coverage_comb = nans(1, num_combination_methods);
        err_rms_binwise_comb = nans(num_combination_methods, num_bins_coarse);
        sigma_rms_binwise_comb = nans(num_combination_methods, num_bins_coarse);
        pdf_w = cell(num_combination_methods, num_uncertainty_methods);

        % ===============================
        %% aggregate weights
        % ===============================
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

        % ===============================
        %% calculate statistics for combined schemes
        % ===============================
        for combination_method_index = 1:num_combination_methods
            %% remove invalid measurements for the combined uncertainty    
            [unc_combined_all_valid{combination_method_index}, valid_indices] = nan_invalid_measurements(real(unc_combined_all{combination_method_index}), min_error_threshold, max_error_threshold);    
            % ===============================
            %% calculate statistics    
            % ===============================    
            % calculate rms of unweighted uncertainty
            sigma_rms_comb(combination_method_index) = rms(unc_combined_all_valid{combination_method_index}, 'omitnan');
            
            %% calculate pdf of the error and uncertainty distributions
            sigma_all_current = real(unc_combined_all_valid{combination_method_index});
            % calculate pdf of unweighted uncertainty
            pdf_sigma_comb(combination_method_index, :) = histcounts(sigma_all_current(isfinite(sigma_all_current)), bins, 'normalization', 'pdf');
            cdf_sigma_comb(combination_method_index, :) = histcounts(sigma_all_current(isfinite(sigma_all_current)), bins, 'normalization', 'cdf');
            
            % ===============================
            %% calculate coverage
            % ===============================
            % unweighted coverage
            coverage_comb(combination_method_index) = calculate_coverage_percentage(err_all_valid, unc_combined_all_valid{combination_method_index});
            
            %% calculate rms error and uncertainty binwise
            
            % calculate binwise rms of error and uncertainty for unweighted
            [err_rms_binwise_comb(combination_method_index, :), sigma_rms_binwise_comb(combination_method_index, :)] = calculate_rms_binwise(err_all_valid, sigma_all_current, max_error_threshold, num_bins_coarse);

            % ===============================
            %% calculate histogram of weights
            % ===============================    
            % calculate pdf of weights for the uncertainty schemes
            for uncertainty_method_index = 1:num_uncertainty_methods
                weights_all_valid{combination_method_index, uncertainty_method_index} = weights_all{combination_method_index, uncertainty_method_index}(valid_indices);

                [pdf_w{combination_method_index, uncertainty_method_index}, ~] = histcounts(weights_all_valid{combination_method_index, uncertainty_method_index}, bins_weights, 'normalization', 'pdf');
                % [pdf_w{combination_method_index, uncertainty_method_index}, ~] = histcounts(weights_all_valid{combination_method_index, uncertainty_method_index}, bins_weights);
            end
                    
        end

        % ===================================
        %% estimate error from uncertainties
        % ===================================
        % seed random number generator
        rng(100);
        % individual methods
        num_methods = num_uncertainty_methods;

        err_est_valid_individual = cell(1, num_methods);
        pdf_err_est_individual = nans(num_methods, numel(bins)-1);
        cdf_err_est_individual = nans(num_methods, numel(bins)-1);

        for method_index = 1:num_methods
            num_valid_trials = numel(sigma_all_valid{method_index});
            err_est_current = randn(1, num_valid_trials) .* sigma_all_valid{method_index};
            % only retain valid indices
            err_est_valid_individual{method_index} = err_est_current(abs(err_est_current) > min_error_threshold & ...
                                        abs(err_est_current) < max_error_threshold);	

            % calculate pdf of estimated error
            [pdf_err_est_individual(method_index, :), ~] = histcounts(abs(err_est_valid_individual{method_index}), bins, 'normalization', 'pdf');
            cdf_err_est_individual(method_index, :) = histcounts(abs(err_est_valid_individual{method_index}), bins, 'normalization', 'cdf');
        end


        % combined methods
        num_methods = num_combination_methods;

        err_est_valid_combined = cell(1, num_methods);
        pdf_err_est_combined = nans(num_methods, numel(bins)-1);
        cdf_err_est_combined = nans(num_methods, numel(bins)-1);

        for method_index = 1:num_methods
            num_valid_trials = numel(unc_combined_all_valid{method_index});
            err_est_current = randn(1, num_valid_trials) .* unc_combined_all_valid{method_index};
            % only retain valid indices
            err_est_valid_combined{method_index} = err_est_current(abs(err_est_current) > min_error_threshold & ...
                                        abs(err_est_current) < max_error_threshold);	

            % calculate pdf of estimated error
            [pdf_err_est_combined(method_index, :), ~] = histcounts(abs(err_est_valid_combined{method_index}), bins, 'normalization', 'pdf');
            cdf_err_est_combined(method_index, :) = histcounts(abs(err_est_valid_combined{method_index}), bins, 'normalization', 'cdf');
        end
        
        % ===============================
        %% calculate histograms for true and estimated error
        % ===============================
        % scaling factor
        scaling_factor = num_trials * (bins(2) - bins(1));

        % error
        h_err = pdf_err * scaling_factor;
        % replace zero values with small finite number
        indices = h_err == 0;
        h_err(indices) = small_value;

        % individual methods
        h_err_est_individual = nans(num_uncertainty_methods, num_bins-1);
        % loop through individual methods
        for uncertainty_method_index = 1:num_uncertainty_methods
            % calcualte histogram from pdf
            h_err_est_individual(uncertainty_method_index, :) = pdf_err_est_individual(uncertainty_method_index, :) * scaling_factor;
            % replace zero values with small finite number
            indices = h_err_est_individual(uncertainty_method_index, :) == 0;
            h_err_est_individual(uncertainty_method_index, indices) = small_value;
        end

        % combined methods
        h_err_est_combined = nans(num_combination_methods, num_bins-1);
        % loop through combination method
        for combination_method_index = 1:num_combination_methods
            % calcualte histogram from pdf
            h_err_est_combined(combination_method_index, :) = pdf_err_est_combined(combination_method_index, :) * scaling_factor;
            % replace zero values with small finite number
            indices = h_err_est_combined(combination_method_index, :) == 0;
            h_err_est_combined(combination_method_index, indices) = small_value;
        end

        % ===============================
        %% calculate histogram distances between true and estimated error
        % ===============================
        % individual methods
        d_err_est_individual = nans(num_uncertainty_methods, num_distance_methods);
        % loop through individual methods
        for uncertainty_method_index = 1:num_uncertainty_methods
            % loop through distance calculation method
            for distance_method_index = 1:num_distance_methods
                % current method
                current_distance_method = histogram_distance_methods{distance_method_index};
                % calculate distance measure
                d_err_est_individual(uncertainty_method_index, distance_method_index) = pdist2(h_err, h_err_est_individual(uncertainty_method_index, :), str2func(current_distance_method));
            end
        end

        % combined methods
        d_err_est_combined = nans(num_combination_methods, num_distance_methods);
        % loop through combination method
        for combination_method_index = 1:num_combination_methods
            % loop through distance calculation method
            for distance_method_index = 1:num_distance_methods
                % current method
                current_distance_method = histogram_distance_methods{distance_method_index};
                % calculate distance measure
                d_err_est_combined(combination_method_index, distance_method_index) = pdist2(h_err, h_err_est_combined(combination_method_index, :), str2func(current_distance_method));
            end
        end
        
        % % ===============================
        % %% display results
        % % ===============================
        for distance_method_index = 1 %:num_distance_methods
            fprintf('==============================\n');
            fprintf('%s\n', histogram_distance_methods{distance_method_index});
            fprintf('Indiv: %.2f, %.2f, %.2f\n', d_err_est_individual(:, distance_method_index));    
            fprintf('Comb: %.2f, %.2f, %.2f\n', d_err_est_combined(:, distance_method_index));    
        end

        % ===============================
        %% aggregate results across all processing
        % ===============================
        % error
        err_all_valid_all_datasets = [err_all_valid_all_datasets, err_all_valid];
        % individual uncertainties
        for method_index = 1:num_uncertainty_methods
            sigma_all_valid_all_datasets{method_index} = [sigma_all_valid_all_datasets{method_index}, sigma_all_valid{method_index}];
            err_est_valid_individual_all_datasets{method_index} = [err_est_valid_individual_all_datasets{method_index}, err_est_valid_individual{method_index}];            
        end
        % combined uncertainties
        for method_index = 1:num_combination_methods
            unc_combined_all_valid_all_datasets{method_index} = [unc_combined_all_valid_all_datasets{method_index}, ...
                                                                 unc_combined_all_valid{method_index}];        
            err_est_valid_combined_all_datasets{method_index} = [err_est_valid_combined_all_datasets{method_index}, ...
                                                                err_est_valid_combined{method_index}];
        end

        % weights
        for combination_method_index = 1:num_combination_methods        
            % aggregate weights
            for uncertainty_method_index = 1:num_uncertainty_methods
                weights_all_valid_all_datasets{combination_method_index, uncertainty_method_index} = [weights_all_valid_all_datasets{combination_method_index, uncertainty_method_index}, ...
                                                                    weights_all_valid{combination_method_index, uncertainty_method_index}];
            end
        end

        % ===================================
        %% save results
        % ===================================
        filename = fullfile(current_write_directory, 'error_uncertainty_statistics.mat');
        save(filename, 'err_all_valid', 'pdf_err', 'cdf_err', 'err_rms', 'coverage_err', ...
            'sigma_all_valid', 'pdf_sigma', 'cdf_sigma', 'sigma_rms', 'coverage', 'err_rms_binwise', 'sigma_rms_binwise', ...
            'unc_combined_all_valid', 'weights_all_valid', 'pdf_sigma_comb', 'cdf_sigma_comb', 'sigma_rms_comb', 'coverage_comb', ...
            'err_rms_binwise_comb', 'sigma_rms_binwise_comb', 'pdf_w', ...
            'h_err', 'h_err_est_individual', 'h_err_est_combined', 'd_err_est_individual', 'd_err_est_combined', ...
            'err_est_valid_individual', 'pdf_err_est_individual', 'cdf_err_est_individual', ...
            'err_est_valid_combined', 'pdf_err_est_combined', 'cdf_err_est_combined');

    end
end

% =============================
%% calculate rms values
% =============================
err_rms_all_datasets = rms(err_all_valid_all_datasets(isfinite(err_all_valid_all_datasets)));
sigma_rms_all_datasets = nans(1, num_uncertainty_methods);
for method_index = 1:num_uncertainty_methods
    sigma_current = sigma_all_valid_all_datasets{method_index};
    sigma_rms_all_datasets(method_index) = rms(sigma_current(isfinite(sigma_current)));
end

sigma_rms_comb_all_datasets = nans(1, num_combination_methods);
for method_index = 1:num_combination_methods
    sigma_current = unc_combined_all_valid_all_datasets{method_index};
    sigma_rms_comb_all_datasets(method_index) = rms(sigma_current(isfinite(sigma_current)));
end

% ==============================
% calculate pdfs
% ==============================
% calculate pdf of error
[pdf_err_all_datasets, ~] = histcounts(err_all_valid_all_datasets(isfinite(err_all_valid_all_datasets)), bins, 'normalization', 'pdf');
cdf_err_all_datasets = histcounts(err_all_valid_all_datasets(isfinite(err_all_valid_all_datasets)), bins, 'normalization', 'cdf');

% calculate pdf of weights
pdf_w_all_datasets = cell(num_combination_methods, num_uncertainty_methods);
for combination_method_index = 1:num_combination_methods        
    for uncertainty_method_index = 1:num_uncertainty_methods
        [pdf_w_all_datasets{combination_method_index, uncertainty_method_index}, ~] = histcounts(weights_all_valid_all_datasets{combination_method_index, uncertainty_method_index}, bins_weights, 'normalization', 'pdf');
    end
end

% ===================================
%% estimate error from uncertainties
% ===================================
% seed random number generator
rng(100);
% individual methods
num_methods = num_uncertainty_methods;
pdf_err_est_individual_all_datasets = nans(num_methods, numel(bins)-1);
cdf_err_est_individual_all_datasets = nans(num_methods, numel(bins)-1);

for method_index = 1:num_methods
    % calculate pdf of estimated error
    [pdf_err_est_individual_all_datasets(method_index, :), ~] = histcounts(abs(err_est_valid_individual_all_datasets{method_index}), bins, 'normalization', 'pdf');
    cdf_err_est_individual_all_datasets(method_index, :) = histcounts(abs(err_est_valid_individual_all_datasets{method_index}), bins, 'normalization', 'cdf');
end

% combined methods
num_methods = num_combination_methods;
pdf_err_est_combined_all_datasets = nans(num_methods, numel(bins)-1);
cdf_err_est_combined_all_datasets = nans(num_methods, numel(bins)-1);

for method_index = 1:num_methods
    % calculate pdf of estimated error
    [pdf_err_est_combined_all_datasets(method_index, :), ~] = histcounts(abs(err_est_valid_combined_all_datasets{method_index}), bins, 'normalization', 'pdf');
    cdf_err_est_combined_all_datasets(method_index, :) = histcounts(abs(err_est_valid_combined_all_datasets{method_index}), bins, 'normalization', 'cdf');
end

% ===============================
%% calculate histograms for true and estimated error
% ===============================
% scaling factor
scaling_factor = num_trials * (bins(2) - bins(1));

% error
h_err_all_datasets = pdf_err_all_datasets * scaling_factor;
% replace zero values with small finite number
indices = h_err_all_datasets == 0;
h_err_all_datasets(indices) = small_value;

% individual methods
h_err_est_individual_all_datasets = nans(num_uncertainty_methods, num_bins-1);
% loop through individual methods
for uncertainty_method_index = 1:num_uncertainty_methods
    % calcualte histogram from pdf
    h_err_est_individual_all_datasets(uncertainty_method_index, :) = pdf_err_est_individual_all_datasets(uncertainty_method_index, :) * scaling_factor;
    % replace zero values with small finite number
    indices = h_err_est_individual_all_datasets(uncertainty_method_index, :) == 0;
    h_err_est_individual_all_datasets(uncertainty_method_index, indices) = small_value;
end

% combined methods
h_err_est_combined_all_datasets = nans(num_combination_methods, num_bins-1);
% loop through combination method
for combination_method_index = 1:num_combination_methods
    % calcualte histogram from pdf
    h_err_est_combined_all_datasets(combination_method_index, :) = pdf_err_est_combined_all_datasets(combination_method_index, :) * scaling_factor;
    % replace zero values with small finite number
    indices = h_err_est_combined_all_datasets(combination_method_index, :) == 0;
    h_err_est_combined_all_datasets(combination_method_index, indices) = small_value;
end

% ===============================
%% calculate histogram distances between true and estimated error
% ===============================
% individual methods
d_err_est_individual_all_datasets = nans(num_uncertainty_methods, num_distance_methods);
% loop through individual methods
for uncertainty_method_index = 1:num_uncertainty_methods
    % loop through distance calculation method
    for distance_method_index = 1:num_distance_methods
        % current method
        current_distance_method = histogram_distance_methods{distance_method_index};
        % calculate distance measure
        d_err_est_individual_all_datasets(uncertainty_method_index, distance_method_index) = pdist2(h_err_all_datasets, h_err_est_individual_all_datasets(uncertainty_method_index, :), str2func(current_distance_method));
    end
end

% combined methods
d_err_est_combined_all_datasets = nans(num_combination_methods, num_distance_methods);
% loop through combination method
for combination_method_index = 1:num_combination_methods
    % loop through distance calculation method
    for distance_method_index = 1:num_distance_methods
        % current method
        current_distance_method = histogram_distance_methods{distance_method_index};
        % calculate distance measure
        d_err_est_combined_all_datasets(combination_method_index, distance_method_index) = pdist2(h_err_all_datasets, h_err_est_combined_all_datasets(combination_method_index, :), str2func(current_distance_method));
    end
end

% ===================================
%% save aggregate results
% ===================================
% directory to store results
current_write_directory = fullfile(top_write_directory, 'all-combined', ...
['trials=' num2str(num_trials, '%d')], ['ppr' num2str(percentage_particles_remove, '%.2f') '_nrs' num2str(num_resampling_trials) '-prewindow'], ...
['max_error=' num2str(max_error_threshold, '%.2f') 'pix-new']);
mkdir_c(current_write_directory);

% save aggregated results
filename = fullfile(current_write_directory, 'error_uncertainty_statistics.mat');
save(filename, 'err_all_valid_all_datasets', 'sigma_all_valid_all_datasets', 'err_est_valid_individual_all_datasets', ...
    'unc_combined_all_valid_all_datasets', 'err_est_valid_combined_all_datasets', ...
    'err_rms_all_datasets', 'sigma_rms_all_datasets', 'sigma_rms_comb_all_datasets', ...
    'pdf_err_all_datasets', 'cdf_err_all_datasets', 'pdf_w_all_datasets', ...
    'd_err_est_individual_all_datasets', 'd_err_est_combined_all_datasets');
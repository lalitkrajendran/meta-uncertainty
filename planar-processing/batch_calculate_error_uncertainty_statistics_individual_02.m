clear
close all
clc

restoredefaultpath;
addpath('../histogram_distance/') 
addpath(genpath('/scratch/shannon/c/aether/Projects/BOS/general-codes/matlab-codes/'));
% addpath('../prana-uncertainty-average-dc-new-im-2-save-deform-cs/');
addpath('../prana/');
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
dataset_name_array = {'PivChal03B'}; % 'PivChal05B'; 'stagnation_flow'; 'Vortex_Ring'; 'Jetdata'};
num_datasets = numel(dataset_name_array);

% array of uncertainty methods
individual_method_array = {'IM'; 'MC'; 'CS'};
num_individual_methods = numel(individual_method_array);

% array of snr metrics
snr_metric_array = {'PPR'; 'MI'};
num_snr_methods = numel(snr_metric_array);
% array of uncertainty combination model to be used
% 'unweighted', 'global-weight-var', 'global-weight-std',
% 'local-weight-var', 'local-weight-std'
% combination_method_array = {'unwt'; 'pd-var'; 'entropy'}; % 'prob'};
combination_method_array = {'unwt'; 'var'; 'entropy'}; % 'prob'};
num_combination_methods = numel(combination_method_array);

% component names
component_names = {'x'; 'y'};
num_components = numel(component_names);

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
% case name
resampling_case_name = ['ppr' num2str(percentage_particles_remove, '%.2f') '_nrs' num2str(num_resampling_trials) '-prewindow-randpart-snr-newproc'];

% ===============================
%% declare variables to store aggregated results
% ===============================
err_all_valid_all_datasets = [];

unc_individual_all_valid_all_datasets = create_empty_cell_array(1, num_individual_methods);
unc_resampled_all_valid_all_datasets = create_empty_cell_array(1, num_individual_methods);
unc_combined_all_valid_all_datasets = create_empty_cell_array(1, num_combination_methods);

snr_individual_all_valid_all_datasets = create_empty_cell_array(1, num_snr_methods);
snr_resampled_all_valid_all_datasets = create_empty_cell_array(1, num_snr_methods);

weights_all_valid_all_datasets = create_empty_cell_array(num_combination_methods, num_individual_methods);

err_est_valid_individual_all_datasets = create_empty_cell_array(1, num_individual_methods);
err_est_valid_combined_all_datasets = create_empty_cell_array(1, num_combination_methods);

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
        ['trials=' num2str(num_trials, '%d')], resampling_case_name);

        % directory to save results for this case
        current_write_directory = fullfile(current_read_directory, ['max_error=' num2str(max_error_threshold, '%.2f') 'pix']);
        mkdir_c(current_write_directory);

        %% load results
        load(fullfile(current_read_directory, 'monte_carlo_results.mat'));
        load(fullfile(current_read_directory, 'combined_uncertainties.mat'));

        % ===============================
        %% aggregate measurements
        % ===============================
        % ---------------
        % error
        % ---------------
        err_all = abs([err_U(:); err_V(:)]);

        % ---------------
        % error - sub
        % ---------------
        err_sub_all = abs([err_U_sub_trials(:); err_V_sub_trials(:)]);

        % ---------------
        % error - resampled
        % ---------------
        err_resampled_all = [];
        for trial_index = 1:num_trials
            err_resampled_all = [err_resampled_all, abs(err_U_resampling{trial_index}), abs(err_V_resampling{trial_index})];
        end

        % ---------------
        % uncertainties - individual
        % ---------------
        unc_individual_all = create_empty_cell_array(1, num_individual_methods);
        for trial_index = 1:num_trials
            if ~isempty(unc_trials{trial_index})
                unc_individual_all{1} = [unc_individual_all{1}, unc_trials{trial_index}.imx, unc_trials{trial_index}.imy];
                unc_individual_all{2} = [unc_individual_all{2}, unc_trials{trial_index}.mcx, unc_trials{trial_index}.mcy];
                unc_individual_all{3} = [unc_individual_all{3}, unc_trials{trial_index}.csx, unc_trials{trial_index}.csy];
            else
                unc_individual_all{1} = [unc_individual_all{1}, NaN, NaN];
                unc_individual_all{2} = [unc_individual_all{2}, NaN, NaN];
                unc_individual_all{3} = [unc_individual_all{3}, NaN, NaN];
            end
        end

        % ---------------
        % uncertainties - individual, calculated on the sub region
        % ---------------
        unc_individual_sub_all = create_empty_cell_array(1, num_individual_methods);
        for trial_index = 1:num_trials
            for method_index = 1:num_individual_methods
                if ~isempty(unc_sub_trials{trial_index})
                    unc_current = [unc_sub_trials{trial_index}{method_index}.x, unc_sub_trials{trial_index}{method_index}.y];
                else
                    unc_current = nans(1, 2);
                end
                unc_individual_sub_all{method_index} = [unc_individual_sub_all{method_index}, unc_current];
            end
        end

        % ---------------
        % uncertainties - resampled
        % ---------------
        unc_resampled_all = {[]; []; []};
        for trial_index = 1:num_trials
            for individual_method_index = 1:num_individual_methods
                unc_current = unc_resampled{trial_index}{individual_method_index};
                if ~isempty(unc_current)
                    unc_resampled_all{individual_method_index} = [unc_resampled_all{individual_method_index}, unc_current.x, unc_current.y];
                end
            end
        end

        % ---------------
        % uncertainties - combined
        % ---------------
        unc_combined_all = {[]; []; []; []};
        for trial_index = 1:num_trials
            unc_combined_current = unc_combined{trial_index};            
            for combination_method_index = 1:num_combination_methods        
                % aggregate uncertainties
                if ~isempty(unc_combined_current)
                    unc_current = [unc_combined_current{combination_method_index}.x, unc_combined_current{combination_method_index}.y];
                else
                    unc_current = [NaN, NaN];
                end                
                unc_combined_all{combination_method_index} = [unc_combined_all{combination_method_index}, unc_current];
            end
        end

        % ---------------
        % weights
        % ---------------
        weights_all = repmat({[]}, 4, 3);
        for trial_index = 1:num_trials
            weights_current = weights{trial_index};            
            for combination_method_index = 1:num_combination_methods                        
                % aggregate weights
                for individual_method_index = 1:num_individual_methods
                    if ~isempty(weights_current)
                        weights_all{combination_method_index, individual_method_index} = [weights_all{combination_method_index, individual_method_index}, weights_current{combination_method_index}.x(individual_method_index), ...
                                                                            weights_current{combination_method_index}.y(individual_method_index)];
                    else
                        weights_all{combination_method_index, individual_method_index} = [weights_all{combination_method_index, individual_method_index}, NaN, NaN];
                    end
                end
            end
        end

        % ---------------
        % snr metric - individual
        % ---------------
        snr_individual_all = create_empty_cell_array(1, num_snr_methods);
        for trial_index = 1:num_trials
            for method_index = 1:num_snr_methods
                if ~isempty(unc_trials{trial_index})
                    snr_current = snr_metric_trials{trial_index}.(snr_metric_array{method_index}) * ones(1, num_components);
                else
                    snr_current = nans(1, num_components);
                end
                snr_individual_all{method_index} = [snr_individual_all{method_index}, snr_current];
            end
        end

        % ---------------
        % snr metric - individual, calculated on the sub region
        % ---------------
        snr_individual_sub_all = create_empty_cell_array(1, num_snr_methods);
        for trial_index = 1:num_trials
            for method_index = 1:num_snr_methods
                if ~isempty(unc_trials{trial_index})
                    snr_current = snr_metric_sub_trials{trial_index}.(snr_metric_array{method_index}) * ones(1, num_components);
                else
                    snr_current = nans(1, num_components);
                end
                snr_individual_sub_all{method_index} = [snr_individual_sub_all{method_index}, snr_current];
            end
        end

        % ---------------
        % snr metric - resampled
        % ---------------
        snr_resampled_all = {[]; []};
        for trial_index = 1:num_trials
            for snr_method_index = 1:num_snr_methods
                snr_current = snr_resampled{trial_index}{snr_method_index};
                if ~isempty(snr_current)
                    snr_resampled_all{snr_method_index} = [snr_resampled_all{snr_method_index}, snr_current];
                end
            end
        end

        % ===============================
        %% identify valid measurements
        % ===============================
        % ---------------
        % error
        % ---------------
        % identify valid indices
        [~, valid_indices_err] = nan_invalid_measurements(err_all, min_error_threshold, max_error_threshold);
        % create logical array with valid indices
        logical_array_err = zeros(1, numel(err_all));
        logical_array_err(valid_indices_err) = 1;

        % ---------------
        % uncertainty - individual
        % ---------------
        logical_array_unc_individual = ones(1, numel(err_all));
        for individual_method_index = 1:num_individual_methods 
            % identify valid indices
            [~, valid_indices_unc_individual{individual_method_index}] = nan_invalid_measurements(real(unc_individual_all{individual_method_index}), min_error_threshold, max_error_threshold);            
            % create logical array with valid indices
            logical_array_temp = zeros(1, numel(unc_individual_all{individual_method_index}));
            logical_array_temp(valid_indices_unc_individual{individual_method_index}) = 1;
    
            % calculate product of all individual logical arrays            
            logical_array_unc_individual = logical_array_unc_individual .* logical_array_temp;
        end

        % ---------------
        % uncertainty - combined
        % ---------------
        logical_array_unc_combined = ones(1, numel(err_all));
        for combination_method_index = 1:num_combination_methods 
            % identify valid indices
            [~, valid_indices_unc_combined{combination_method_index}] = nan_invalid_measurements(real(unc_combined_all{combination_method_index}), min_error_threshold, max_error_threshold);    
            % create logical array with valid indices
            logical_array_temp = zeros(1, numel(unc_combined_all{combination_method_index}));
            logical_array_temp(valid_indices_unc_combined{combination_method_index}) = 1;    

            % calculate product of all combined logical arrays
            logical_array_unc_combined = logical_array_unc_combined .* logical_array_temp;
        end

        % array of common valid indices
        logical_array_common = logical_array_err .* logical_array_unc_individual .* logical_array_unc_combined;
        % identify valid common indices
        valid_indices_common = find(logical_array_common);
        % number of valid trials
        num_valid_trials = numel(valid_indices_common);

        % ===============================
        %% separate valid measurements
        % ===============================
        % ---------------
        % error
        % ---------------
        err_all_valid = err_all(valid_indices_common);

        % ---------------
        % error - resampled
        % ---------------
        err_resampled_all_valid = [];
        err_resampled_all_valid = extract_valid_resampled_elements(err_resampled_all, valid_indices_common, num_resampling_trials);

        % ---------------
        % uncertainty - individual
        % ---------------
        unc_individual_all_valid = cell(1, num_individual_methods);
        for individual_method_index = 1:num_individual_methods 
            unc_individual_all_valid{individual_method_index} = unc_individual_all{individual_method_index}(valid_indices_common);
        end

        % ---------------
        % uncertainty - resampled
        % ---------------
        unc_resampled_all_valid = cell(1, num_individual_methods);
        for individual_method_index = 1:num_individual_methods
            unc_resampled_all_valid{individual_method_index} = extract_valid_resampled_elements(unc_resampled_all{individual_method_index}, valid_indices_common, num_resampling_trials);
        end

        % ---------------
        % uncertainty - combined
        % ---------------
        unc_combined_all_valid = cell(1, num_combination_methods);
        for combination_method_index = 1:num_combination_methods 
            unc_combined_all_valid{combination_method_index} = unc_combined_all{combination_method_index}(valid_indices_common);
        end

        % ---------------
        % weights
        % ---------------
        weights_all_valid = repmat({[]}, 4, 3);
        for combination_method_index = 1:num_combination_methods
            for individual_method_index = 1:num_individual_methods
                weights_all_valid{combination_method_index, individual_method_index} = weights_all{combination_method_index, individual_method_index}(valid_indices_common);
            end
        end
        
        % ---------------
        % snr - individual
        % ---------------
        snr_individual_all_valid = cell(1, num_snr_methods);
        for snr_method_index = 1:num_snr_methods
            snr_individual_all_valid{snr_method_index} = snr_individual_all{snr_method_index}(valid_indices_common);
        end

        % ---------------
        % snr - resampled
        % ---------------
        snr_resampled_all_valid = cell(1, num_snr_methods);
        for snr_method_index = 1:num_snr_methods
            snr_resampled_all_valid{snr_method_index} = extract_valid_resampled_elements(snr_resampled_all{snr_method_index}, valid_indices_common, num_resampling_trials);
        end        

        % ===============================
        %% calculate resampling means
        % ===============================
        % ---------------
        % uncertainty - resampled
        % ---------------
        mean_unc_resampled = nans(num_individual_methods, num_valid_trials);
        for individual_method_index = 1:num_individual_methods 
            mean_unc_resampled(individual_method_index, :) = calculate_resampled_means(unc_resampled_all_valid{individual_method_index}, num_resampling_trials);
        end

        % ---------------
        % snr metric - resampled
        % ---------------
        mean_snr_resampled = nans(num_snr_methods, num_valid_trials);
        for snr_method_index = 1:num_snr_methods
            mean_snr_resampled(snr_method_index, :) = calculate_resampled_means(snr_resampled_all_valid{snr_method_index}, num_resampling_trials); 
        end
        
        % ===============================
        %% calculate rms
        % ===============================
        % ---------------
        % error
        % ---------------
        rms_err = rms(err_all_valid);

        % ---------------
        % error - resampled
        % ---------------
        rms_err_resampled = rms(err_resampled_all_valid);

        % ---------------
        % uncertainty - individual
        % ---------------
        rms_unc_individual = nans(1, num_individual_methods);
        for individual_method_index = 1:num_individual_methods 
            rms_unc_individual(individual_method_index) = rms(unc_individual_all_valid{individual_method_index}, 'omitnan');
        end

        % ---------------
        % uncertainty - resampled
        % ---------------
        rms_unc_resampled = nans(1, num_individual_methods);
        for individual_method_index = 1:num_individual_methods 
            rms_unc_resampled(individual_method_index) = rms(unc_resampled_all_valid{individual_method_index}, 'omitnan');
        end

        % ---------------
        % uncertainty - combined
        % ---------------
        rms_unc_combined = nans(1, num_combination_methods);
        for combination_method_index = 1:num_combination_methods 
            rms_unc_combined(combination_method_index) = rms(unc_combined_all_valid{combination_method_index});
        end

        % ===============================
        %% calculate pdf, cdf
        % ===============================
        % ---------------
        % error
        % ---------------
        pdf_err = histcounts(err_all_valid, bins, 'normalization', 'pdf');
        cdf_err = histcounts(err_all_valid, bins, 'normalization', 'cdf');

        % ---------------
        % error - resampled
        % ---------------
        pdf_err_resampled = histcounts(err_resampled_all_valid, bins, 'normalization', 'pdf');
        cdf_err_resampled = histcounts(err_resampled_all_valid, bins, 'normalization', 'cdf');

        % allocate memory
        pdf_err_resampled_trials = nans(num_trials, num_bins-1);
        % calculate pdf of resampled errors for each trial
        for trial_index = 1:num_valid_trials
            % index of the first element for this trial
            start_index = (trial_index - 1) * num_resampling_trials + 1;
            % index of the last element for this trial
            stop_index = trial_index * num_resampling_trials;
            % extract resampled values for this trial
            err_current_trial = err_resampled_all_valid(start_index:stop_index);
            % calculate pdf of errors uncertainties
            pdf_err_resampled_trials(trial_index, :) = histcounts(err_current_trial(isfinite(err_current_trial)), bins, 'normalization', 'pdf');                
        end

        % ---------------
        % uncertainty - individual
        % ---------------
        pdf_unc_individual = nans(num_individual_methods, num_bins-1);
        cdf_unc_individual = nans(num_individual_methods, num_bins-1);
        for individual_method_index = 1:num_individual_methods 
            % extract uncertainties
            unc_current = real(unc_individual_all_valid{individual_method_index});            
            % only retain finite values
            unc_current = unc_current(isfinite(unc_current));
            % calculate pdf
            pdf_unc_individual(individual_method_index, :) = histcounts(unc_current, bins, 'normalization', 'pdf');
            % calculate cdf
            cdf_unc_individual(individual_method_index, :) = histcounts(unc_current, bins, 'normalization', 'cdf');
        end

        % ---------------
        % uncertainty - resampled
        % ---------------
        pdf_unc_resampled = nans(num_individual_methods, num_bins-1);
        pdf_unc_resampled_trials = cell(1, num_individual_methods);
        cdf_unc_resampled = nans(num_individual_methods, num_bins-1);
        for individual_method_index = 1:num_individual_methods 
            % extract uncertainties
            unc_current = real(unc_resampled_all_valid{individual_method_index});            
            % only retain finite values
            unc_current = unc_current(isfinite(unc_current));
            % calculate pdf
            pdf_unc_resampled(individual_method_index, :) = histcounts(unc_current, bins, 'normalization', 'pdf');
            % calculate cdf
            cdf_unc_resampled(individual_method_index, :) = histcounts(unc_current, bins, 'normalization', 'cdf');

            % allocate memory
            pdf_unc_resampled_trials{individual_method_index} = nans(num_trials, num_bins-1);
            % calculate pdf of resampled uncertainties for each trial
            for trial_index = 1:numel(valid_indices_common)
                % index of the first element for this trial
                start_index = (trial_index - 1) * num_resampling_trials + 1;
                % index of the last element for this trial
                stop_index = trial_index * num_resampling_trials;
                % extract resampled uncertainties for this trial
                unc_current_trial = unc_resampled_all_valid{individual_method_index}(start_index:stop_index);
                % calculate pdf of resampled uncertainties
                pdf_unc_resampled_trials{individual_method_index}(trial_index, :) = histcounts(unc_current_trial(isfinite(unc_current_trial)), bins, 'normalization', 'pdf');                
            end
        end

        % ---------------
        % uncertainty - combined
        % ---------------
        pdf_unc_combined = nans(num_combination_methods, num_bins-1);
        cdf_unc_combined = nans(num_combination_methods, num_bins-1);
        for combination_method_index = 1:num_combination_methods
            % extract uncertainties
            unc_current = real(unc_combined_all_valid{combination_method_index});
            % only retain finite values
            unc_current = unc_current(isfinite(unc_current));
            % calculate pdf
            pdf_unc_combined(combination_method_index, :) = histcounts(unc_current(isfinite(unc_current)), bins, 'normalization', 'pdf');
            % calculate cdf
            cdf_unc_combined(combination_method_index, :) = histcounts(unc_current(isfinite(unc_current)), bins, 'normalization', 'cdf');                            
        end

        % ---------------
        % weights
        % ---------------
        pdf_w = cell(num_combination_methods, num_individual_methods);
        for combination_method_index = 1:num_combination_methods
            for individual_method_index = 1:num_individual_methods
                [pdf_w{combination_method_index, individual_method_index}, ~] = histcounts(weights_all_valid{combination_method_index, individual_method_index}, bins_weights, 'normalization', 'pdf');
            end
        end

        % ===================================
        %% estimate error from uncertainties
        % ===================================
        % seed random number generator
        rng(100);

        % ---------------
        % individual 
        % ---------------
        num_methods = num_individual_methods;

        err_est_valid_individual = cell(1, num_methods);
        pdf_err_est_individual = nans(num_methods, numel(bins)-1);
        cdf_err_est_individual = nans(num_methods, numel(bins)-1);

        for method_index = 1:num_methods
            num_valid_trials = numel(unc_individual_all_valid{method_index});
            err_est_current = randn(1, num_valid_trials) .* unc_individual_all_valid{method_index};
            % only retain valid indices
            err_est_valid_individual{method_index} = err_est_current(abs(err_est_current) > min_error_threshold & ...
                                        abs(err_est_current) < max_error_threshold);	

            % calculate pdf of estimated error
            [pdf_err_est_individual(method_index, :), ~] = histcounts(abs(err_est_valid_individual{method_index}), bins, 'normalization', 'pdf');
            cdf_err_est_individual(method_index, :) = histcounts(abs(err_est_valid_individual{method_index}), bins, 'normalization', 'cdf');
        end

        % ---------------
        % combined 
        % ---------------
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

        % ---------------
        % true error
        % ---------------
        h_err = pdf_err * scaling_factor;
        % replace zero values with small finite number
        indices = h_err == 0;
        h_err(indices) = small_value;

        % ---------------
        % estimated error - individual 
        % ---------------
        h_err_est_individual = nans(num_individual_methods, num_bins-1);
        % loop through individual methods
        for individual_method_index = 1:num_individual_methods
            % calcualte histogram from pdf
            h_err_est_individual(individual_method_index, :) = pdf_err_est_individual(individual_method_index, :) * scaling_factor;
            % replace zero values with small finite number
            indices = h_err_est_individual(individual_method_index, :) == 0;
            h_err_est_individual(individual_method_index, indices) = small_value;
        end

        % ---------------
        % estimated error - combined
        % ---------------
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
        % ---------------
        % individual 
        % ---------------
        d_err_est_individual = nans(num_individual_methods, num_distance_methods);
        % loop through individual methods
        for individual_method_index = 1:num_individual_methods
            % loop through distance calculation method
            for distance_method_index = 1:num_distance_methods
                % current method
                current_distance_method = histogram_distance_methods{distance_method_index};
                % calculate distance measure
                d_err_est_individual(individual_method_index, distance_method_index) = pdist2(h_err, h_err_est_individual(individual_method_index, :), str2func(current_distance_method));
            end
        end

        % ---------------
        % combined 
        % ---------------
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

        % ===================================
        %% save results
        % ===================================
        filename = fullfile(current_write_directory, 'error_uncertainty_statistics_02.mat');
        save(filename, 'unc_individual_all', 'unc_individual_sub_all', 'snr_individual_sub_all', ...
            'valid_indices_common', 'err_all_valid', 'pdf_err', 'cdf_err', 'rms_err', ...
            'err_resampled_all_valid', 'pdf_err_resampled', 'cdf_err_resampled', 'rms_err_resampled', 'pdf_err_resampled_trials', ...
            'unc_individual_all_valid', 'pdf_unc_individual', 'cdf_unc_individual', 'rms_unc_individual', ...
            'unc_resampled_all_valid', 'pdf_unc_resampled', 'cdf_unc_resampled', 'rms_unc_resampled', 'pdf_unc_resampled_trials', ...
            'unc_combined_all_valid', 'weights_all_valid', 'pdf_unc_combined', 'cdf_unc_combined', 'rms_unc_combined', ...            
            'weights_all_valid', 'pdf_w', ...
            'snr_individual_all_valid', 'snr_resampled_all_valid', 'mean_unc_resampled', 'mean_snr_resampled',...
            'h_err', 'h_err_est_individual', 'h_err_est_combined', 'd_err_est_individual', 'd_err_est_combined', ...
            'err_est_valid_individual', 'pdf_err_est_individual', 'cdf_err_est_individual', ...
            'err_est_valid_combined', 'pdf_err_est_combined', 'cdf_err_est_combined');

        return;
        % ===============================
        %% aggregate results across all processing
        % ===============================
        % ---------------
        % error
        % ---------------
        err_all_valid_all_datasets = [err_all_valid_all_datasets, err_all_valid'];

        % ---------------
        % uncertainties - individual
        % ---------------
        for method_index = 1:num_individual_methods
            unc_individual_all_valid_all_datasets{method_index} = [unc_individual_all_valid_all_datasets{method_index}, unc_individual_all_valid{method_index}];
            err_est_valid_individual_all_datasets{method_index} = [err_est_valid_individual_all_datasets{method_index}, err_est_valid_individual{method_index}];            
        end

        % ---------------
        % uncertainties - resampled
        % ---------------
        rms_unc_resampled = nans(1, num_individual_methods);
        for method_index = 1:num_individual_methods 
            unc_resampled_all_valid_all_datasets{method_index} = [unc_resampled_all_valid_all_datasets{method_index}, unc_resampled_all_valid{method_index}];
        end

        % ---------------
        % uncertainties - combined
        % ---------------
        for method_index = 1:num_combination_methods
            unc_combined_all_valid_all_datasets{method_index} = [unc_combined_all_valid_all_datasets{method_index}, unc_combined_all_valid{method_index}];        
            err_est_valid_combined_all_datasets{method_index} = [err_est_valid_combined_all_datasets{method_index}, err_est_valid_combined{method_index}];
        end

        % ----------
        % weights
        % ----------
        for combination_method_index = 1:num_combination_methods        
            for individual_method_index = 1:num_individual_methods
                weights_all_valid_all_datasets{combination_method_index, individual_method_index} = [weights_all_valid_all_datasets{combination_method_index, individual_method_index}, ...
                                                                    weights_all_valid{combination_method_index, individual_method_index}];
            end
        end

        % ---------------
        % snr - individual
        % ---------------
        for method_index = 1:num_snr_methods
            snr_individual_all_valid_all_datasets{method_index} = [snr_individual_all_valid_all_datasets{method_index}, snr_individual_all_valid{method_index}];
        end

        % ---------------
        % snr - resampled
        % ---------------
        for method_index = 1:num_snr_methods 
            snr_resampled_all_valid_all_datasets{method_index} = [snr_resampled_all_valid_all_datasets{method_index}, snr_resampled_all_valid{method_index}];
        end
    end    
end

% =============================
%% calculate rms values
% =============================
% ---------------
% error
% ---------------
rms_err_all_datasets = rms(err_all_valid_all_datasets(isfinite(err_all_valid_all_datasets)), 'omitnan');

% ---------------
% uncertainty - individual
% ---------------
rms_unc_individual_all_datasets = nans(1, num_individual_methods);
for method_index = 1:num_individual_methods
    sigma_current = unc_individual_all_valid_all_datasets{method_index};
    rms_unc_individual_all_datasets(method_index) = rms(sigma_current(isfinite(sigma_current)), 'omitnan');
end

% ---------------
% uncertainty - resampled
% ---------------
rms_unc_resampled_all_datasets = nans(1, num_individual_methods);
for individual_method_index = 1:num_individual_methods 
    rms_unc_resampled_all_datasets(individual_method_index) = rms(unc_resampled_all_valid_all_datasets{individual_method_index}, 'omitnan');
end

% ---------------
% uncertainty - combined
% ---------------
rms_unc_combined_all_datasets = nans(1, num_combination_methods);
for method_index = 1:num_combination_methods
    sigma_current = unc_combined_all_valid_all_datasets{method_index};
    rms_unc_combined_all_datasets(method_index) = rms(sigma_current(isfinite(sigma_current)), 'omitnan');
end

% ==============================
% calculate pdf and cdf
% ==============================
% ---------------
% error
% ---------------
pdf_err_all_datasets = histcounts(err_all_valid_all_datasets(isfinite(err_all_valid_all_datasets)), bins, 'normalization', 'pdf');
cdf_err_all_datasets = histcounts(err_all_valid_all_datasets(isfinite(err_all_valid_all_datasets)), bins, 'normalization', 'cdf');

% ---------------
% uncertainty - individual
% ---------------
pdf_unc_individual_all_datasets = nans(num_individual_methods, num_bins-1);
cdf_unc_individual_all_datasets = nans(num_individual_methods, num_bins-1);
for individual_method_index = 1:num_individual_methods 
    % extract uncertainties
    unc_current = real(unc_individual_all_valid_all_datasets{individual_method_index});            
    % only retain finite values
    unc_current = unc_current(isfinite(unc_current));
    % calculate pdf
    pdf_unc_individual_all_datasets(individual_method_index, :) = histcounts(unc_current, bins, 'normalization', 'pdf');
    % calculate cdf
    cdf_unc_individual_all_datasets(individual_method_index, :) = histcounts(unc_current, bins, 'normalization', 'cdf');
end

% ---------------
% uncertainty - resampled
% ---------------
pdf_unc_resampled_all_datasets = nans(num_individual_methods, num_bins-1);
cdf_unc_resampled_all_datasets = nans(num_individual_methods, num_bins-1);
for individual_method_index = 1:num_individual_methods 
    % extract uncertainties
    unc_current = real(unc_resampled_all_valid_all_datasets{individual_method_index});            
    % only retain finite values
    unc_current = unc_current(isfinite(unc_current));
    % calculate pdf
    pdf_unc_resampled_all_datasets(individual_method_index, :) = histcounts(unc_current, bins, 'normalization', 'pdf');
    % calculate cdf
    cdf_unc_resampled_all_datasets(individual_method_index, :) = histcounts(unc_current, bins, 'normalization', 'cdf');
end

% ---------------
% uncertainty - combined
% ---------------
pdf_unc_combined_all_datasets = nans(num_combination_methods, num_bins-1);
cdf_unc_combined_all_datasets = nans(num_combination_methods, num_bins-1);
for combination_method_index = 1:num_combination_methods
    % extract uncertainties
    unc_current = real(unc_combined_all_valid_all_datasets{combination_method_index});
    % only retain finite values
    unc_current = unc_current(isfinite(unc_current));
    % calculate pdf
    pdf_unc_combined_all_datasets(combination_method_index, :) = histcounts(unc_current(isfinite(unc_current)), bins, 'normalization', 'pdf');
    % calculate cdf
    cdf_unc_combined_all_datasets(combination_method_index, :) = histcounts(unc_current(isfinite(unc_current)), bins, 'normalization', 'cdf');                            
end

% ---------------
% weights
% ---------------
pdf_w_all_datasets = cell(num_combination_methods, num_individual_methods);
for combination_method_index = 1:num_combination_methods        
    for individual_method_index = 1:num_individual_methods
        [pdf_w_all_datasets{combination_method_index, individual_method_index}, ~] = histcounts(weights_all_valid_all_datasets{combination_method_index, individual_method_index}, bins_weights, 'normalization', 'pdf');
    end
end

% ===================================
%% estimate error from uncertainties
% ===================================
% seed random number generator
rng(100);

% ---------------
% individual methods
% ---------------
num_methods = num_individual_methods;
pdf_err_est_individual_all_datasets = nans(num_methods, numel(bins)-1);
cdf_err_est_individual_all_datasets = nans(num_methods, numel(bins)-1);

for method_index = 1:num_methods
    % calculate pdf of estimated error
    [pdf_err_est_individual_all_datasets(method_index, :), ~] = histcounts(abs(err_est_valid_individual_all_datasets{method_index}), bins, 'normalization', 'pdf');
    cdf_err_est_individual_all_datasets(method_index, :) = histcounts(abs(err_est_valid_individual_all_datasets{method_index}), bins, 'normalization', 'cdf');
end

% ---------------
% combined methods
% ---------------
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

% ---------------
% error
% ---------------
h_err_all_datasets = pdf_err_all_datasets * scaling_factor;
% replace zero values with small finite number
indices = h_err_all_datasets == 0;
h_err_all_datasets(indices) = small_value;

% ---------------
% individual methods
% ---------------
h_err_est_individual_all_datasets = nans(num_individual_methods, num_bins-1);
% loop through individual methods
for individual_method_index = 1:num_individual_methods
    % calcualte histogram from pdf
    h_err_est_individual_all_datasets(individual_method_index, :) = pdf_err_est_individual_all_datasets(individual_method_index, :) * scaling_factor;
    % replace zero values with small finite number
    indices = h_err_est_individual_all_datasets(individual_method_index, :) == 0;
    h_err_est_individual_all_datasets(individual_method_index, indices) = small_value;
end

% ---------------
% combined methods
% ---------------
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
% ---------------
% individual methods
% ---------------
d_err_est_individual_all_datasets = nans(num_individual_methods, num_distance_methods);
% loop through individual methods
for individual_method_index = 1:num_individual_methods
    % loop through distance calculation method
    for distance_method_index = 1:num_distance_methods
        % current method
        current_distance_method = histogram_distance_methods{distance_method_index};
        % calculate distance measure
        d_err_est_individual_all_datasets(individual_method_index, distance_method_index) = pdist2(h_err_all_datasets, h_err_est_individual_all_datasets(individual_method_index, :), str2func(current_distance_method));
    end
end

% ---------------
% combined methods
% ---------------
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
['trials=' num2str(num_trials, '%d')], resampling_case_name, ...
['max_error=' num2str(max_error_threshold, '%.2f') 'pix-new']);
mkdir_c(current_write_directory);

% save aggregated results
filename = fullfile(current_write_directory, 'error_uncertainty_statistics_02.mat');
save(filename, 'err_all_valid_all_datasets', 'unc_individual_all_valid_all_datasets', 'unc_resampled_all_valid_all_datasets', 'unc_combined_all_valid_all_datasets', ... 
                'rms_err_all_datasets', 'rms_unc_individual_all_datasets', 'rms_unc_resampled_all_datasets', 'rms_unc_combined_all_datasets', ...
                'pdf_err_all_datasets', 'cdf_err_all_datasets', 'pdf_unc_individual_all_datasets', 'cdf_unc_individual_all_datasets', ... 
                'pdf_unc_resampled_all_datasets', 'cdf_unc_resampled_all_datasets', 'pdf_unc_combined_all_datasets', 'cdf_unc_combined_all_datasets', ...
                'pdf_w_all_datasets', 'snr_individual_all_valid_all_datasets', 'snr_resampled_all_valid_all_datasets', ...
                'err_est_valid_individual_all_datasets', 'err_est_valid_combined_all_datasets', ...    
                'd_err_est_individual_all_datasets', 'd_err_est_combined_all_datasets');
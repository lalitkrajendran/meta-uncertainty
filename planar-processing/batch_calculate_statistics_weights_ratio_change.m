% This script performs resampling calculations for a range of particle removal
% percentage and plots the change in uncertainties for different methods

clear
close all
clc

restoredefaultpath;
addpath(genpath('/scratch/shannon/c/aether/Projects/BOS/general-codes/matlab-codes/'));
addpath('../general-codes/')
addpath('../histogram_distance/')
setup_default_settings;

% ============================
%% read/write settings
% ============================
% window resolution
window_resolution_array = [64, 32];
num_window_resolution = numel(window_resolution_array);
% pass number
pass_number = 4;
% directory containing images
top_image_directory = fullfile('/scratch/shannon/a/aether/Projects/PIV_PTV_Uncertainty_Quantification/Planar_uncertainty_work/Images/');
% directory containing files to be read
top_read_directory = fullfile('/scratch/shannon/a/aether/Projects/PIV_PTV_Uncertainty_Quantification/uncertainty-combination/analysis/results/planar-dataset', 'experiment-new');
% array of data set names
dataset_name_array = {'PivChal03B'; 'PivChal05B'; 'stagnation_flow'; 'Vortex_Ring'; 'Jetdata'};
num_datasets = numel(dataset_name_array);
% array of uncertainty methods
individual_method_array = {'IM'; 'MC'; 'CS'};
num_individual_methods = numel(individual_method_array);
% array of snr metrics
snr_metric_array = {'PPR'; 'MI'};
num_snr_methods = numel(snr_metric_array);
% components
components = {'x'; 'y'};
num_components = numel(components);
% directory where results of this analysis are to be saved
top_write_directory = fullfile('/scratch/shannon/a/aether/Projects/PIV_PTV_Uncertainty_Quantification/uncertainty-combination/analysis/results/planar-dataset', 'monte-carlo', 'individual-datasets');
mkdir_c(top_write_directory);
% number of trials
num_trials = 1e3;

% ============================
%% resampling settings
% ============================
% resampling method names
resampling_method_names = {'remove-paired'; 'remove-random'; 'add-random'};
% resampling_method_names = {'add-random'};
num_resampling_methods = numel(resampling_method_names);
% number of particles to remove
percentage_particles_remove_array = 0.05:0.05:0.25;
num_ppr = numel(percentage_particles_remove_array);
% number of resampling trials
num_resampling_trials = 1e2; %3;
% display_resampled_images? (true/false)
display_resampled_images = 0;
% method names for plotting
resampling_method_names_plot = {'Removing Paired Particles'; 'Removing Random Particles'; 'Adding Random Particles'};
resampling_method_names_plot_short = {'R-P', 'R-R', 'A-R'};
% % case name
% resampling_case_name = ['ppr' num2str(percentage_particles_remove, '%.2f') '_nrs' num2str(num_resampling_trials) '-prewindow-randpart-snr-newproc'];

% ============================
%% particle identification settings
% ============================
% intensity threshold for particle identification
intensity_threshold = 100;
% particle diameter (pix.)
d_p = 3;
% particle sizing settings
sizeprops = struct;
sizeprops.method = 'IWC';
sizeprops.p_area = 2; % this will be modified for each image
sizeprops.sigma = 4;
sizeprops.errors = 0;
sizeprops.save_dir = [];
% display id results?
display_id_results = 0;

% ========================
%% weight settings
% ========================
% number of bins for weights
num_bins_w = 30;
% bins for weights
bins_w = linspace(0, 1, num_bins_w);

% ========================
%% statistics settings
% ========================
min_error_threshold = 0;
max_error_threshold = 1; %00; %0.2;
num_bins = 30;
bins = linspace(-max_error_threshold, max_error_threshold, num_bins*2);
bins_abs = linspace(1e-3, max_error_threshold, num_bins);
% methods to calculate histogram distances
histogram_distance_method = 'total_variation_distance';
% categories
histogram_distance_categories = categorical({individual_method_array{:}, resampling_method_names_plot_short{:}});
% sort in the desired order
histogram_distance_categories = reordercats(histogram_distance_categories, ...
                                {individual_method_array{:}, resampling_method_names_plot_short{:}});
% error model ('gaussian' or 'lognormal')
% error_models = {'gaussian'; 'lognormal'; 'exp'};
% error_models = {'rayleigh'; 'lognormal'; 'exp'};
error_models = {'gaussian'};
num_models = numel(error_models);

% metric type ('rms', 'iqr', 'std')
metric_name = 'iqr';

% % ============================
% %% pre-load all dataset
% % ============================
% fprintf('Loading all datasets into memory\n');

% results_all = cell(num_window_resolution, num_datasets);
% jobfile_all = cell(num_window_resolution, num_datasets);
% files_im1 = cell(num_window_resolution, num_datasets);
% files_im2 = cell(num_window_resolution, num_datasets);

% % loop through window resolutions
% for window_resolution_index = 1:num_window_resolution
%     % loop through datasets
%     for dataset_index = 1:num_datasets
%         % name of the data set
%         dataset_name = dataset_name_array{dataset_index};        
%         fprintf('Dataset: %s\n', dataset_name);

%         % ============================
%         %% load data
%         % ============================
%         % results directory for current data set
%         current_results_directory = fullfile(top_read_directory, dataset_name, ['WS' num2str(window_resolution_index)]);

%         % directory containing vectors
%         vectors_directory = fullfile(current_results_directory, 'vectors-new');

%         % load results for vectors, errors and uncertainties
%         results_all{window_resolution_index, dataset_index} = load_directory_data(vectors_directory, ['*pass' num2str(pass_number) '*.mat']);

%         % load jobfile
%         jobfile_all{window_resolution_index, dataset_index} = load(fullfile(current_results_directory, 'jobfile.mat'));
        
%         % ============================
%         %% load listing of deformed images
%         % ============================        
%         % directory containing deformed images for current data set
%         deformed_images_directory = fullfile(vectors_directory, 'imDeform');
        
%         % get list of im1 files in the directory
%         files_im1{window_resolution_index, dataset_index} = get_directory_listing(deformed_images_directory, 'PIV*im1d*.mat');
%         % get list of im2 files in the directory
%         files_im2{window_resolution_index, dataset_index} = get_directory_listing(deformed_images_directory, 'PIV*im2d*.mat');        
%     end
% end

% % ============================
% %% load errors for all datasets
% % ============================
% fprintf('Loading all errors into memory\n');
% errors_all = cell(num_window_resolution, num_datasets);

% % loop through window resolutions
% for window_resolution_index = 1:num_window_resolution
%     % loop through datasets
%     for dataset_index = 1:num_datasets
%         % name of the data set        
%         dataset_name = dataset_name_array{dataset_index};
        
%         fprintf('Dataset: %s\n', dataset_name);
        
%         %% Load data

%         % results directory for current data set
%         current_results_directory = fullfile(top_read_directory, dataset_name, ['WS' num2str(window_resolution_index)]);

%         % load results for vectors, errors and uncertainties
%         errors_all{window_resolution_index, dataset_index} = load(fullfile(current_results_directory, 'errors-new.mat'));        
%     end
% end

% initialize variables
err_all_consolidated = [];
unc_indiv_all_consolidated = [];
unc_comb_all_consolidated = [];
valid_trials_consolidated = [];
num_valid_trials_datasets = [];
err_est_indiv_consolidated = [];
err_est_comb_consolidated = [];
wt_consolidated = cell(1, num_resampling_methods);

% ============================
%% loop through window resolutions
% ============================
for window_resolution_index = 1:num_window_resolution
    fprintf('WS: %d\n', window_resolution_index);

    % ============================
    %% loop through datasets
    % ============================
    for dataset_index = 1:num_datasets
        % dataset name
        dataset_name = dataset_name_array{dataset_index};
        fprintf('dataset: %s\n', dataset_name);

        % ================================================
        % load results
        % ================================================
        % directory to save results for this case
        current_read_directory = fullfile(top_write_directory, dataset_name, ['WS' num2str(window_resolution_index)], 'weights-rms-change-study');
        % file name
        filename = ['resampling-statistics-nt=' num2str(num_trials, '%d') 'nr='  num2str(num_resampling_trials, '%d') '.mat'];
        % load results to file
        load(fullfile(current_read_directory, filename));
        % file name for weights and combined uncertainties
        filename = ['resampling-statistics-nt=' num2str(num_trials, '%d') 'nr='  num2str(num_resampling_trials, '%d') '-weights-unc-' metric_name '-02-new.mat'];
        % load results to file
        load(fullfile(current_read_directory, filename));

        % ==========================
        % calculate error/uncertainty statistics
        % ==========================
        % extract valid measurements
        [valid_trials, err_valid, unc_indiv_valid, unc_comb_valid] = extract_valid_errors_uncertainties(err, unc_indiv, unc_comb, components, ...
                                                                        min_error_threshold * ones(1, num_components), max_error_threshold * ones(1, num_components));
        num_valid_trials = sum(double(valid_trials));

        % display valid trials
        fprintf('%--------------------------------\n');
        fprintf('num_valid_trials: %d, num_trials: %d\n', num_valid_trials, num_trials);
        err_temp = abs([err.x; err.y]);
        unc_indiv_temp = abs([unc_indiv.x(:); unc_indiv.y(:)]);
        unc_comb_temp = abs([unc_comb.x(:); unc_comb.y(:)]);
        fprintf('Error: %f, %f\n', min(err_temp), max(err_temp));
        fprintf('Unc, Indiv: %f, %f\n', min(unc_indiv_temp), max(unc_indiv_temp));
        fprintf('Unc, Comb: %f, %f\n',  min(unc_comb_temp), max(unc_comb_temp));        
        fprintf('%--------------------------------\n');

        % aggregate measurements across components
        [err_all, unc_indiv_all, unc_comb_all] = aggregate_errors_uncertainties(err_valid, unc_indiv_valid, unc_comb_valid, components);

        % calculate rms
        [err_rms, unc_indiv_rms, unc_comb_rms, err_avg, confint_err] = calculate_error_uncertainty_rms(err_all, unc_indiv_all, unc_comb_all);

        % calculate pdf of weights
        pdf_w = calculate_pdf_weights(wt_trials, bins_w, valid_trials, num_individual_methods, num_resampling_methods, components);

        % estimate error from uncertainties
        [err_est_indiv, err_est_comb] = estimate_error_from_uncertainty_gaussian(err_all, unc_indiv_all, unc_comb_all);

        % calculate histogram distance between true and estimated error
        [d_err_est_indiv, d_err_est_comb] = calculate_histogram_distance_err_true_est(err_all, err_est_indiv, err_est_comb, bins, histogram_distance_method);

        % % ==========================
        % % save results
        % % ==========================
        % % directory to save results for this case
        % current_read_directory = fullfile(top_write_directory, dataset_name, ['WS' num2str(window_resolution_index)], 'weights-rms-change-study');
        % % file name
        % filename = ['resampling-statistics-nt=' num2str(num_trials, '%d') 'nr='  num2str(num_resampling_trials, '%d') '-statistics-' metric_name '-new-thresh=' num2str(max_error_threshold, '%.2f') 'pix.mat'];
        % % load results to file
        % save(fullfile(current_read_directory, filename), 'valid_trials', 'num_valid_trials', 'pdf_w', ...
        %                                                 'unc_indiv_all', 'unc_indiv_rms', 'unc_comb_all', 'unc_comb_rms', ...
        %                                                 'confint_err', 'err_all', 'err_rms', ...
        %                                                 'err_est_indiv', 'err_est_comb', ...
        %                                                 'd_err_est_indiv', 'd_err_est_comb');

        % ==========================
        % consolidate results across datasets
        % ==========================
        err_all_consolidated = [err_all_consolidated; err_all];
        unc_indiv_all_consolidated = [unc_indiv_all_consolidated; unc_indiv_all];
        unc_comb_all_consolidated = [unc_comb_all_consolidated; unc_comb_all];
        valid_trials_consolidated = [valid_trials_consolidated; valid_trials];
        num_valid_trials_datasets = [num_valid_trials_datasets; num_valid_trials];

        for resampling_method_index = 1:num_resampling_methods
            if dataset_index == 1 && window_resolution_index == 1
                wt_consolidated{resampling_method_index}.x = wt_trials{resampling_method_index}.x;
                wt_consolidated{resampling_method_index}.y = wt_trials{resampling_method_index}.y;
            else
                wt_consolidated{resampling_method_index}.x = [wt_consolidated{resampling_method_index}.x; wt_trials{resampling_method_index}.x];
                wt_consolidated{resampling_method_index}.y = [wt_consolidated{resampling_method_index}.y; wt_trials{resampling_method_index}.y];                    
            end    
        end

        err_est_indiv_consolidated = [err_est_indiv_consolidated; err_est_indiv];
        err_est_comb_consolidated = [err_est_comb_consolidated; err_est_comb];
    end
end

% ==========================
% calculate consolidated statistics
% ==========================

% number of valid trials overall
num_valid_trials_consolidated = sum(double(valid_trials_consolidated));
valid_trials_consolidated = logical(valid_trials_consolidated);

% calculate rms
[err_rms_consolidated, unc_indiv_rms_consolidated, unc_comb_rms_consolidated, ~, confint_err_consolidated] = calculate_error_uncertainty_rms(err_all_consolidated', unc_indiv_all_consolidated, unc_comb_all_consolidated);

% pdf of weights
pdf_w_consolidated = calculate_pdf_weights(wt_consolidated, bins_w, valid_trials_consolidated, num_individual_methods, num_resampling_methods, components);

% calculate histogram distance between true and estimated error
[d_err_est_indiv_consolidated, d_err_est_comb_consolidated] = calculate_histogram_distance_err_true_est(err_all_consolidated, err_est_indiv_consolidated, err_est_comb_consolidated, ...
                                                    bins, histogram_distance_method);
return;
% ==========================
% load global properties
% ==========================
% directory to save results for this case
write_directory = fullfile(top_write_directory, 'weights-rms-change-study-consolidated');
% file name
filename = ['resampling-statistics-nt=' num2str(num_trials, '%d') 'nr='  num2str(num_resampling_trials, '%d') '-physical-prop.mat'];
% load results
load(fullfile(write_directory, filename));

% ==========================
% calculate uncertainty as function of fractional displacement
% ==========================
num_bins_disp = 10;
bins_disp = linspace(-0.5, 0.5, num_bins_disp);
frac_U = U_valid_consolidated - round(U_valid_consolidated);
[indices_disp, num_indices_disp, err_rms_binned_disp, unc_indiv_rms_binned_disp, unc_comb_rms_binned_disp] = calculate_binwise_rms(frac_U, err_all_consolidated, unc_indiv_all_consolidated, ...
                                                                                    unc_comb_all_consolidated, bins_disp);
frac_U_ref = U_ref_valid_consolidated - round(U_ref_valid_consolidated);
[indices_disp_ref, num_indices_disp_ref, err_rms_binned_disp_ref, unc_indiv_rms_binned_disp_ref, unc_comb_rms_binned_disp_ref] = calculate_binwise_rms(frac_U_ref, err_all_consolidated, unc_indiv_all_consolidated, ...
                                                                                    unc_comb_all_consolidated, bins_disp);

% ==========================
% calculate uncertainty as function of shear
% ==========================
% number of bins to be used
num_bins_grad = 10;
% bin edges
bins_grad = linspace(0, 0.025, num_bins_grad);

grad_U = abs(U_shear_valid_consolidated);
[indices_grad, num_indices_grad, err_rms_binned_grad, unc_indiv_rms_binned_grad, unc_comb_rms_binned_grad] = calculate_binwise_rms(grad_U, err_all_consolidated, unc_indiv_all_consolidated, ...
                                                                                    unc_comb_all_consolidated, bins_grad);
grad_U_ref = abs(U_shear_ref_valid_consolidated);
[indices_grad_ref, num_indices_grad_ref, err_rms_binned_grad_ref, unc_indiv_rms_binned_grad_ref, unc_comb_rms_binned_grad_ref] = calculate_binwise_rms(grad_U_ref, err_all_consolidated, unc_indiv_all_consolidated, ...
                                                                                    unc_comb_all_consolidated, bins_grad);

% ==========================
% save consolidated results
% ==========================
% directory to save results for this case
write_directory = fullfile(top_write_directory, 'weights-rms-change-study-consolidated');
mkdir_c(write_directory);
% file name
filename = ['resampling-statistics-nt=' num2str(num_trials, '%d') 'nr='  num2str(num_resampling_trials, '%d') '-statistics-' metric_name '-new-thresh=' num2str(max_error_threshold, '%.2f') 'pix.mat'];
% load results to file
save(fullfile(write_directory, filename), 'valid_trials_consolidated', 'num_valid_trials_consolidated', 'num_valid_trials_datasets', ...
                                                'err_all_consolidated', 'unc_indiv_all_consolidated', 'unc_comb_all_consolidated', ...
                                                'pdf_w_consolidated', 'err_est_indiv_consolidated', 'err_est_comb_consolidated', ...
                                                'err_rms_consolidated', 'confint_err_consolidated', ...
                                                'unc_indiv_rms_consolidated', 'unc_comb_rms_consolidated',...                                                
                                                'err_est_indiv_consolidated', 'err_est_comb_consolidated', ...
                                                'd_err_est_indiv_consolidated', 'd_err_est_comb_consolidated', ...
                                                'bins_disp', 'bins_grad', ...
                                                'frac_U', 'indices_disp', 'num_indices_disp', 'err_rms_binned_disp', 'unc_indiv_rms_binned_disp', 'unc_comb_rms_binned_disp', ...
                                                'frac_U_ref', 'indices_disp_ref', 'num_indices_disp_ref', 'err_rms_binned_disp_ref', 'unc_indiv_rms_binned_disp_ref', 'unc_comb_rms_binned_disp_ref', ...
                                                'grad_U', 'indices_grad', 'num_indices_grad', 'err_rms_binned_grad', 'unc_indiv_rms_binned_grad', 'unc_comb_rms_binned_grad', ...
                                                'grad_U_ref', 'indices_grad_ref', 'num_indices_grad_ref', 'err_rms_binned_grad_ref', 'unc_indiv_rms_binned_grad_ref', 'unc_comb_rms_binned_grad_ref');


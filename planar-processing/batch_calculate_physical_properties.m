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

% ============================
%% pre-load all dataset
% ============================
fprintf('Loading all datasets into memory\n');

results_all = cell(num_window_resolution, num_datasets);
jobfile_all = cell(num_window_resolution, num_datasets);
files_im1 = cell(num_window_resolution, num_datasets);
files_im2 = cell(num_window_resolution, num_datasets);

% loop through window resolutions
for window_resolution_index = 1:num_window_resolution
    % loop through datasets
    for dataset_index = 1:num_datasets
        % name of the data set
        dataset_name = dataset_name_array{dataset_index};        
        fprintf('Dataset: %s\n', dataset_name);

        % ============================
        %% load data
        % ============================
        % results directory for current data set
        current_results_directory = fullfile(top_read_directory, dataset_name, ['WS' num2str(window_resolution_index)]);

        % directory containing vectors
        vectors_directory = fullfile(current_results_directory, 'vectors-new');

        % load results for vectors, errors and uncertainties
        results_all{window_resolution_index, dataset_index} = load_directory_data(vectors_directory, ['*pass' num2str(pass_number) '*.mat']);

        % load jobfile
        jobfile_all{window_resolution_index, dataset_index} = load(fullfile(current_results_directory, 'jobfile.mat'));
        
        % ============================
        %% load listing of deformed images
        % ============================        
        % directory containing deformed images for current data set
        deformed_images_directory = fullfile(vectors_directory, 'imDeform');
        
        % get list of im1 files in the directory
        files_im1{window_resolution_index, dataset_index} = get_directory_listing(deformed_images_directory, 'PIV*im1d*.mat');
        % get list of im2 files in the directory
        files_im2{window_resolution_index, dataset_index} = get_directory_listing(deformed_images_directory, 'PIV*im2d*.mat');        
    end
end

% ============================
%% load errors for all datasets
% ============================
fprintf('Loading all errors into memory\n');
errors_all = cell(num_window_resolution, num_datasets);

% loop through window resolutions
for window_resolution_index = 1:num_window_resolution
    % loop through datasets
    for dataset_index = 1:num_datasets
        % name of the data set        
        dataset_name = dataset_name_array{dataset_index};
        
        fprintf('Dataset: %s\n', dataset_name);
        
        %% Load data

        % results directory for current data set
        current_results_directory = fullfile(top_read_directory, dataset_name, ['WS' num2str(window_resolution_index)]);

        % load results for vectors, errors and uncertainties
        errors_all{window_resolution_index, dataset_index} = load(fullfile(current_results_directory, 'errors-new.mat'));        
    end
end

% --------------------------
% initialize variables
% --------------------------
U_consolidated = [];
U_ref_consolidated = [];

U_valid_consolidated = [];
U_ref_valid_consolidated = [];

U_grad_consolidated = [];
U_grad_ref_consolidated = [];

U_grad_valid_consolidated = [];
U_grad_ref_valid_consolidated = [];

U_shear_consolidated = [];
U_shear_ref_consolidated = [];

U_shear_valid_consolidated = [];
U_shear_ref_valid_consolidated = [];

U_strain_consolidated = [];
U_strain_ref_consolidated = [];

U_strain_valid_consolidated = [];
U_strain_ref_valid_consolidated = [];

err_rms_consolidated = [];
err_rms_valid_consolidated = [];

err_random_consolidated = [];
err_random_valid_consolidated = [];

fprintf('Calculating global properties\n');

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
        % file name for statistics
        filename = ['resampling-statistics-nt=' num2str(num_trials, '%d') 'nr='  num2str(num_resampling_trials, '%d') '-statistics-' metric_name '-new-thresh=' num2str(max_error_threshold, '%.2f') 'pix.mat'];
        % load results to file
        load(fullfile(current_read_directory, filename));

        % ================================================
        % calculate displacement gradients
        % ================================================
        results = results_all{window_resolution_index, dataset_index};
        errors = errors_all{window_resolution_index, dataset_index};

        % --------------------------
        % initialize
        % --------------------------
        U_trials = nans(num_trials, 1);
        V_trials = nans(num_trials, 1);
        U_ref_trials = nans(num_trials, 1);
        V_ref_trials = nans(num_trials, 1);
        
        U_grad_trials = nans(num_trials, 1);
        V_grad_trials = nans(num_trials, 1);
        U_grad_ref_trials = nans(num_trials, 1);
        V_grad_ref_trials = nans(num_trials, 1);
        
        U_shear_trials = nans(num_trials, 1);
        V_shear_trials = nans(num_trials, 1);
        U_shear_ref_trials = nans(num_trials, 1);
        V_shear_ref_trials = nans(num_trials, 1);

        U_strain_trials = nans(num_trials, 1);
        V_strain_trials = nans(num_trials, 1);
        U_strain_ref_trials = nans(num_trials, 1);
        V_strain_ref_trials = nans(num_trials, 1);

        err_rms_U = nans(num_trials, 1);
        err_rms_V = nans(num_trials, 1);
        
        err_random_U = nans(num_trials, 1);
        err_random_V = nans(num_trials, 1);
        
        parfor trial_index = 1:num_trials

            % extract array indices
            r = r_trials(trial_index);
            c = c_trials(trial_index);

            % --------------------------
            % processed result
            % --------------------------
            % extract displacement value at the current grid point
            U_trials(trial_index) = results{snapshot_trials(trial_index)}.U(r, c);
            V_trials(trial_index) = results{snapshot_trials(trial_index)}.V(r, c);

            % calculate displacement gradient
            [dU_dX, dU_dY, dV_dX, dV_dY] = calculate_displacement_gradient_2D(results{snapshot_trials(trial_index)}.X, results{snapshot_trials(trial_index)}.Y, ...
                                                                            results{snapshot_trials(trial_index)}.U, results{snapshot_trials(trial_index)}.V);

            % calculate absolute magnitude of the gradient
            U_grad = abs(dU_dX) + abs(dU_dY);
            V_grad = abs(dV_dX) + abs(dV_dY);

            % extract value at the current grid point
            U_grad_trials(trial_index) = U_grad(r, c);
            V_grad_trials(trial_index) = V_grad(r, c);

            % calculate shear
            U_shear_trials(trial_index) = dU_dY(r, c);
            V_shear_trials(trial_index) = dV_dX(r, c);

            % calculate strain
            U_strain_trials(trial_index) = dU_dX(r, c);
            V_strain_trials(trial_index) = dV_dY(r, c);

            % extract true displacement value at the current grid point
            U_ref_trials(trial_index) = errors.Ut_interp(r, c, snapshot_trials(trial_index));
            V_ref_trials(trial_index) = errors.Vt_interp(r, c, snapshot_trials(trial_index));

            % --------------------------
            % true solution
            % --------------------------
            % calculate true displacement gradient
            [dU_dX_ref, dU_dY_ref, dV_dX_ref, dV_dY_ref] = calculate_displacement_gradient_2D(errors.X, errors.Y, errors.Ut_interp(:, :, snapshot_trials(trial_index)), ...
                                                                                            errors.Vt_interp(:, :, snapshot_trials(trial_index)));
            
            % calculate absolute magnitude of the gradient
            U_grad_ref = abs(dU_dX_ref) + abs(dU_dY_ref);
            V_grad_ref = abs(dV_dX_ref) + abs(dV_dY_ref);

            % extract true value at the current grid point
            U_grad_ref_trials(trial_index) = U_grad_ref(r, c);
            V_grad_ref_trials(trial_index) = V_grad_ref(r, c);

            % calculate shear
            U_shear_ref_trials(trial_index) = dU_dY_ref(r, c);
            V_shear_ref_trials(trial_index) = dV_dX_ref(r, c);

            % calculate strain
            U_strain_ref_trials(trial_index) = dU_dX_ref(r, c);
            V_strain_ref_trials(trial_index) = dV_dY_ref(r, c);

            % calculate rms of error at this point across snapshots
            err_rms_U(trial_index) = errors.total_U(r, c);
            err_rms_V(trial_index) = errors.total_V(r, c);

            % calculate random error at this point across snapshots
            err_random_U(trial_index) = errors.random_U(r, c);
            err_random_V(trial_index) = errors.random_V(r, c);
        end

        % ==========================
        % extract valid results
        % ==========================        
        U_trials_valid = U_trials(valid_trials);
        V_trials_valid = V_trials(valid_trials);

        U_ref_trials_valid = U_ref_trials(valid_trials);
        V_ref_trials_valid = V_ref_trials(valid_trials);

        U_grad_trials_valid = U_grad_trials(valid_trials);
        V_grad_trials_valid = V_grad_trials(valid_trials);
        U_grad_ref_trials_valid = U_grad_ref_trials(valid_trials);
        V_grad_ref_trials_valid = V_grad_ref_trials(valid_trials);
        
        U_shear_trials_valid = U_shear_trials(valid_trials);
        V_shear_trials_valid = V_shear_trials(valid_trials);
        U_shear_ref_trials_valid = U_shear_ref_trials(valid_trials);
        V_shear_ref_trials_valid = V_shear_ref_trials(valid_trials);

        U_strain_trials_valid = U_strain_trials(valid_trials);
        V_strain_trials_valid = V_strain_trials(valid_trials);
        U_strain_ref_trials_valid = U_strain_ref_trials(valid_trials);
        V_strain_ref_trials_valid = V_strain_ref_trials(valid_trials);

        err_rms_U_valid = err_rms_U(valid_trials);
        err_rms_V_valid = err_rms_V(valid_trials);
        err_random_U_valid = err_random_U(valid_trials);
        err_random_V_valid = err_random_V(valid_trials);
        
        % ==========================
        % save results
        % ==========================
        % directory to save results for this case
        current_read_directory = fullfile(top_write_directory, dataset_name, ['WS' num2str(window_resolution_index)], 'weights-rms-change-study');
        % file name
        filename = ['resampling-statistics-nt=' num2str(num_trials, '%d') 'nr='  num2str(num_resampling_trials, '%d') '-physical-prop.mat'];
        % load results to file
        save(fullfile(current_read_directory, filename), 'U_trials', 'V_trials', 'U_grad_trials', 'V_grad_trials', ...
                                                        'U_shear_trials', 'V_shear_trials', 'U_strain_trials', 'V_strain_trials', ...
                                                        'U_ref_trials', 'V_ref_trials', 'U_grad_ref_trials', 'V_grad_ref_trials', ...
                                                        'U_shear_ref_trials', 'V_shear_ref_trials', 'U_strain_ref_trials', 'V_strain_ref_trials', ...
                                                        'U_trials_valid', 'V_trials_valid', 'U_grad_trials_valid', 'V_grad_trials_valid', ...
                                                        'U_shear_trials_valid', 'V_shear_trials_valid', 'U_strain_trials_valid', 'V_strain_trials_valid', ...
                                                        'U_ref_trials_valid', 'V_ref_trials_valid', 'U_grad_ref_trials_valid', 'V_grad_ref_trials_valid', ...
                                                        'U_shear_ref_trials_valid', 'V_shear_ref_trials_valid', 'U_strain_ref_trials_valid', 'V_strain_ref_trials_valid', ...
                                                        'err_rms_U', 'err_rms_V', 'err_rms_U_valid', 'err_rms_V_valid', ...
                                                        'err_random_U', 'err_random_V', 'err_random_U_valid', 'err_random_V_valid');

        % ==========================
        % consolidate results across datasets
        % ==========================
        U_consolidated = [U_consolidated; U_trials; V_trials];
        U_ref_consolidated = [U_ref_consolidated; U_ref_trials; V_ref_trials];
        
        U_grad_consolidated = [U_grad_consolidated; U_grad_trials; V_grad_trials];
        U_grad_ref_consolidated = [U_grad_ref_consolidated; U_grad_ref_trials; V_grad_ref_trials];
        
        U_shear_consolidated = [U_shear_consolidated; U_shear_trials; V_shear_trials];
        U_strain_consolidated = [U_strain_consolidated; U_strain_trials; V_strain_trials];
        
        U_shear_ref_consolidated = [U_shear_ref_consolidated; U_shear_ref_trials; V_shear_ref_trials];
        U_strain_ref_consolidated = [U_strain_ref_consolidated; U_strain_ref_trials; V_strain_ref_trials];
        
        U_valid_consolidated = [U_valid_consolidated; U_trials_valid; V_trials_valid];
        U_ref_valid_consolidated = [U_ref_valid_consolidated; U_ref_trials_valid; V_ref_trials_valid];

        U_grad_valid_consolidated = [U_grad_valid_consolidated; U_grad_trials_valid; V_grad_trials_valid];
        U_grad_ref_valid_consolidated = [U_grad_ref_valid_consolidated; U_grad_ref_trials_valid; V_grad_ref_trials_valid];

        U_shear_valid_consolidated = [U_shear_valid_consolidated; U_shear_trials_valid; V_shear_trials_valid];
        U_strain_valid_consolidated = [U_strain_valid_consolidated; U_strain_trials_valid; V_strain_trials_valid];
        
        U_shear_ref_valid_consolidated = [U_shear_ref_valid_consolidated; U_shear_ref_trials_valid; V_shear_ref_trials_valid];
        U_strain_ref_valid_consolidated = [U_strain_ref_valid_consolidated; U_strain_ref_trials_valid; V_strain_ref_trials_valid];

        % err_rms_consolidated = [err_rms_consolidated; err_rms_U; err_rms_V];
        % err_rms_valid_consolidated = [err_rms_valid_consolidated; err_rms_U_valid; err_rms_V_valid];

        % err_random_consolidated = [err_random_consolidated; err_random_U; err_random_V];
        % err_random_valid_consolidated = [err_random_valid_consolidated; err_random_U_valid; err_random_V_valid];

    end
end

% ==========================
% save consolidated results
% ==========================
% directory to save results for this case
write_directory = fullfile(top_write_directory, 'weights-rms-change-study-consolidated');
mkdir_c(write_directory);
% file name
filename = ['resampling-statistics-nt=' num2str(num_trials, '%d') 'nr='  num2str(num_resampling_trials, '%d') '-physical-prop.mat'];
% load results to file
save(fullfile(write_directory, filename), 'U_consolidated', 'U_ref_consolidated', 'U_grad_consolidated', 'U_grad_ref_consolidated', ...
                                            'U_valid_consolidated', 'U_ref_valid_consolidated', 'U_grad_valid_consolidated', 'U_grad_ref_valid_consolidated', ...
                                            'U_shear_consolidated', 'U_shear_ref_consolidated', 'U_shear_valid_consolidated', 'U_shear_ref_valid_consolidated', ...
                                            'U_strain_consolidated', 'U_strain_ref_consolidated', 'U_strain_valid_consolidated', 'U_strain_ref_valid_consolidated'); %, ...
                                            % 'err_rms_consolidated', 'err_rms_valid_consolidated', ...
                                            % 'err_random_consolidated', 'err_random_valid_consolidated');

shut_down_parpool();
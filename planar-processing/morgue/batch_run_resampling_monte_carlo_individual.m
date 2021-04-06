% This script calculates the pdf of uncertainties estimated by three
% different algorithms for experimental images by jack-knife resampling.

%%
clear
close all
clc

restoredefaultpath;
addpath(genpath('/scratch/shannon/c/aether/Projects/BOS/general-codes/matlab-codes/'));
% addpath('../prana-uncertainty-average-dc-new-im-2-save-deform-cs/');
addpath ../prana/
addpath ../general-codes/
setup_default_settings;

% ============================
%% read/write settings
% ============================
% window resolution
window_resolution_array = [32, 64];
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
num_uncertainty_methods = numel(individual_method_array);
% directory where results of this analysis are to be saved
top_write_directory = fullfile('/scratch/shannon/a/aether/Projects/PIV_PTV_Uncertainty_Quantification/uncertainty-combination/analysis/results/planar-dataset', 'monte-carlo', 'individual-datasets');
mkdir_c(top_write_directory);
% number of trials
num_trials = 1e2;

% ============================
%% resampling settings
% ============================
% number of particles to remove
% num_particles_remove_array = [5, 10, 15];
percentage_particles_remove = 0.25;
% number of resampling trials
num_resampling_trials = 1e3;
% display_resampled_images? (true/false)
display_resampled_images = 0;
% case name
resampling_case_name = ['ppr' num2str(percentage_particles_remove, '%.2f') '_nrs' num2str(num_resampling_trials) '-prewindow-randpart-snr-newproc'];

% ============================
%% particle identification settings
% ============================
% intensity threshold for particle identification
intensity_threshold = 10;
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

% ============================
%% pre-load all dataset
% ============================
fprintf('Loading all datasets into memory\n');

results_all = cell(num_window_resolution, num_datasets);
jobfile_all = cell(num_window_resolution, num_datasets);
files_im1 = cell(num_window_resolution, num_datasets);
files_im2 = cell(num_window_resolution, num_datasets);

% loop through window resolutions
for window_size_index = 1:num_window_resolution
    % loop through datasets
    for dataset_index = 1:numel(dataset_name_array)
        % name of the data set
        dataset_name = dataset_name_array{dataset_index};        
        fprintf('Dataset: %s\n', dataset_name);

        % ============================
        %% load data
        % ============================
        % results directory for current data set
        current_results_directory = fullfile(top_read_directory, dataset_name, ['WS' num2str(window_size_index)]);

        % directory containing vectors
        vectors_directory = fullfile(current_results_directory, 'vectors-new');

        % load results for vectors, errors and uncertainties
        results_all{window_size_index, dataset_index} = load_directory_data(vectors_directory, ['*pass' num2str(pass_number) '*.mat']);

        % load jobfile
        jobfile_all{window_size_index, dataset_index} = load(fullfile(current_results_directory, 'jobfile.mat'));
        
        % ============================
        %% load listing of deformed images
        % ============================        
        % directory containing deformed images for current data set
        deformed_images_directory = fullfile(vectors_directory, 'imDeform');
        
        % get list of im1 files in the directory
        files_im1{window_size_index, dataset_index} = get_directory_listing(deformed_images_directory, 'PIV*im1d*.mat');
        % get list of im2 files in the directory
        files_im2{window_size_index, dataset_index} = get_directory_listing(deformed_images_directory, 'PIV*im2d*.mat');        
    end
end

% ============================
%% load errors for all datasets
% ============================
fprintf('Loading all errors into memory\n');
errors_all = cell(num_window_resolution, num_datasets);

% loop through window resolutions
for window_size_index = 1:num_window_resolution
    % loop through datasets
    for dataset_index = 1:numel(dataset_name_array)
        % name of the data set        
        dataset_name = dataset_name_array{dataset_index};
        
        fprintf('Dataset: %s\n', dataset_name);
        
        %% Load data

        % results directory for current data set
        current_results_directory = fullfile(top_read_directory, dataset_name, ['WS' num2str(window_size_index)]);

        % load results for vectors, errors and uncertainties
        errors_all{window_size_index, dataset_index} = load(fullfile(current_results_directory, 'errors-new.mat'));        
    end
end

fprintf('running monte-carlo\n');
% set seed for random number generator
rng(0);

% start timer
tic
% ============================
%% loop through window resolutions
% ============================
for window_size_index = 1:num_window_resolution
    fprintf('WS: %d\n', window_size_index);
    % ============================
    %% loop through datasets
    % ============================
    for dataset_index = 1:num_datasets
        % dataset name
        dataset_name = dataset_name_array{dataset_index};
        fprintf('dataset: %s\n', dataset_name);
        % directory to save results for this case
        current_write_directory = fullfile(top_write_directory, dataset_name, ['WS' num2str(window_size_index)], ...
        ['trials=' num2str(num_trials, '%d')], resampling_case_name);
        mkdir_c(current_write_directory);

        % load errors
        % errors = errors_all{window_size_index, dataset_index};
        errors = errors_all{window_size_index, dataset_index};

        % calculate random number for data set
        % dataset_index = randi(num_datasets);
        dataset_index = dataset_index;

        % calculate number of snapshots available for the current dataset and
        % window size
        num_snapshots = numel(results_all{window_size_index, dataset_index});
        if num_snapshots > size(errors.err_U, 3)
            num_snapshots = size(errors.err_U, 3);
        end

        % ====================================
        % extract job file properties
        % ====================================
        % extract job file
        jobfile = jobfile_all{window_size_index, dataset_index}.Data;
        % extract pass settings
        current_pass_settings = jobfile.(['PIV' num2str(pass_number)]);
        % extract window size
        [window_size_x, window_size_y] = extract_window_size(current_pass_settings);
        % extract window resolution
        [window_resolution_x, window_resolution_y] = extract_window_resolution(current_pass_settings);        
        % extract grid resolution
        [grid_resolution_x, grid_resolution_y] = extract_grid_resolution(current_pass_settings);

        % ================================================
        %% determine indices to be used for this dataset
        % ================================================
        if strcmpi(dataset_name, 'Vortex_Ring')        
            % extract co-ordinate arrays for the processing grid
            xvec = errors.X(1, :);
            yvec = errors.Y(:, 1);
            
            % find vector locations in the processing results for which
            % true solution is available
            c = find(xvec > errors.X(1, 1) & xvec < errors.X(1, end));
            r = find(yvec > errors.Y(1, 1) & yvec < errors.Y(end, 1));
        elseif strcmp(dataset_name, 'Jetdata')
            % index limits for processing
            % X 320 to 392 (72); Y 205 to 325 (120)
            cmin_p = 5;
            cmax_p = 23; %303+18-1=320  303+90-1=392  X direction
            rmin_p = 5;
            rmax_p = 35; %188+18-1=205  188+138-1=325 Y direction
                        
            % only retain true solution in the index limits for
            % processing
            c = find(errors.X(1, :) >= 320  & errors.X(1, :) <= 392);
            r = find(errors.Y(:, 1) >= 205  & errors.Y(:, 1) <= 325);
        else
            c = 1:size(errors.X, 2);
            r = 1:size(errors.X, 1);
        end

        % ================================================
        % remove edge points
        % ================================================

        % calculate number of grid points for have a window resolution
        grid_point_buffer_x = 0.5 * window_resolution_x/grid_resolution_x;
        grid_point_buffer_y = 0.5 * window_resolution_y/grid_resolution_y;

        % update row and column list
        c = c(grid_point_buffer_x:(end-grid_point_buffer_x));
        r = r(grid_point_buffer_y:(end-grid_point_buffer_y));

        % ================================================
        %% allocate memory
        % ================================================
        fprintf('allocating memory for variables\n');
        % dataset address
        dataset_index_array = nans(1, num_trials);
        window_size_index_array = nans(1, num_trials);

        % errors
        err_U = nans(1, num_trials);
        err_V = nans(1, num_trials);
        err_U_sub_trials = nans(1, num_trials);
        err_V_sub_trials = nans(1, num_trials);
        err_U_resampling = cell(1, num_trials);
        err_V_resampling = cell(1, num_trials);

        % uncertainties
        unc_trials = cell(1, num_trials);
        unc_sub_trials = cell(1, num_trials);

        % snr metric
        snr_metric_trials = cell(1, num_trials);
        snr_metric_sub_trials = cell(1, num_trials);

        % uncertainty from resampling
        unc_resampling_avg = cell(1, num_trials);
        unc_resampling_std = cell(1, num_trials);
        unc_resampling = cell(1, num_trials);
        snr_metric_resampling = cell(1, num_trials);

        % ================================================
        %% loop through trials and perform analysis
        % ================================================
        for trial_index = 1:num_trials
            fprintf('trial_index: %d\n', trial_index);
            % record dataset index and window resolution index
            window_size_index_array(trial_index) = window_size_index;
            dataset_index_array(trial_index) = dataset_index;

            % calculate random number for the snapshot
            snapshot_index = randi(num_snapshots);

            % load dataset
            results = results_all{window_size_index, dataset_index}{snapshot_index};

            % calculate size of the results array
            [num_rows, num_cols] = size(results.X);
            
            % select a random grid point
            r_trial = randsample(r, 1);
            c_trial = randsample(c, 1);

            % ================================================
            %% extract processing results for this grid point
            % ================================================
            % extract co-ordinates
            X = results.X(r_trial, c_trial);
            Y = results.Y(r_trial, c_trial);
            
            % extract individual uncertainties for this grid point
            unc_current = extract_planar_uncertainties(results.uncertainty2D, r_trial, c_trial);

            % extract snr metric
            snr_metric_current = extract_snr_metric(results.SNRmetric, r_trial, c_trial);

            % ================================================
            %% extract errors for this grid point
            % ================================================
            
            % identify indices on the true solution grid corresponding to the
            % current grid point
            if strcmp(dataset_name, 'Vortex_Ring')
                r_t = find(errors.Y(:, 1) == results.Y(r_trial, 1));
                c_t = find(errors.X(1, :) == results.X(1, c_trial));
            elseif strcmp(dataset_name, 'Jetdata')
                % indices to be used in the true solution
                % rt = 8:2:69;
                % ct = 7:2:43;
                % r_t = r_trial - rmin_p + 1; %find(errors.Y(:, 1) == results.Y(r_trial, 1) + 187);
                % c_t = c_trial - cmin_p + 1; %find(errors.X(1, :) == results.X(1, c_trial) + 302);
                r_t = r_trial;
                c_t = c_trial;
            else
                r_t = r_trial;
                c_t = c_trial;
            end
            
            % extract errors
            err_U_current = errors.err_U(r_t, c_t, snapshot_index);
            err_V_current = errors.err_V(r_t, c_t, snapshot_index);

            % ================================================
            % load deformed image
            % ================================================
            % im1 = imread(fullfile(files_im1{window_size_index, dataset_index}(snapshot_index).folder, files_im1{window_size_index, dataset_index}(snapshot_index).name));
            % im1d = double(flipud(im1));
            saved_img = load(fullfile(files_im1{window_size_index, dataset_index}(snapshot_index).folder, files_im1{window_size_index, dataset_index}(snapshot_index).name));
            im1d = saved_img.im1d;
                
            % im2 = imread(fullfile(files_im2{window_size_index, dataset_index}(snapshot_index).folder, files_im2{window_size_index, dataset_index}(snapshot_index).name));
            % im2d = double(flipud(im2));
            saved_img = load(fullfile(files_im2{window_size_index, dataset_index}(snapshot_index).folder, files_im2{window_size_index, dataset_index}(snapshot_index).name));        
            im2d = saved_img.im2d;
                
            % extract image height and width
            if trial_index == 1
                [image_height, image_width] = size(im1d);
            end

            % extract interrogation window
            % im1_sub = extract_interrogation_window(results, im1d, [window_size_x, window_size_y], r_trial, c_trial);
            % im2_sub = extract_interrogation_window(results, im2d, [window_size_x, window_size_y], r_trial, c_trial);
            im1_sub = extract_interrogation_window(im1d, X, Y, [window_size_x, window_size_y]);
            im2_sub = extract_interrogation_window(im2d, X, Y, [window_size_x, window_size_y]);
            
            % calculate planar uncertainty on the original interrogation windows
            % [unc_sub_current, snr_metric_sub_current] = calculate_planar_uncertainties(im1_sub, im2_sub, [window_size_x, window_size_y], [window_resolution_x, window_resolution_y], individual_method_array, uncertainty_flags);
            % [unc_sub_current, snr_metric_sub_current] = calculate_planar_uncertainties(im1_sub, im2_sub, [window_size_x, window_size_y], [window_resolution_x, window_resolution_y], individual_method_array, uncertainty_flags);
            % [unc_sub_current, snr_metric_sub_current, err_sub_current] = calculate_planar_uncertainties(im1_sub, im2_sub, current_pass_settings, individual_method_array);
            [unc_sub_current, snr_metric_sub_current, U_sub_current, V_sub_current] = calculate_planar_uncertainties(im1_sub, im2_sub, current_pass_settings, individual_method_array);
            err_U_sub_current = U_sub_current + results.Ubs(r_trial, c_trial) - errors.Ut_interp(r_t, c_t, snapshot_index);
            err_V_sub_current = V_sub_current + results.Vbs(r_trial, c_trial) - errors.Vt_interp(r_t, c_t, snapshot_index);

            % ================================================
            %% remove particles and perform resampling calc
            % ================================================            
            % calculate effective particle diameter
            d_p = results.uncertainty2D.Autod(r_trial, c_trial)/sqrt(2);

            % identify particles
            % [x1, y1, x2, y2] = identify_particles(im1_sub, im2_sub, d_p, sizeprops, display_id_results);
            [x1, y1] = identify_particles_02(im1_sub, d_p, sizeprops, display_id_results);
            [x2, y2] = identify_particles_02(im2_sub, d_p, sizeprops, display_id_results);

            % calculate number of particles identified
            num_particles = numel(x1);
            % calculate number of particles to be removed
            num_particles_remove = round(percentage_particles_remove * num_particles);
            % calculate uncertainty from resampling
            % [unc_resampling{trial_index}, snr_metric_resampling{trial_index}] = calculate_uncertainty_resampling(im1_sub, im2_sub, x1, y1, x2, y2, ...
            %                     d_p*ones(size(x1)), num_particles_remove, num_resampling_trials, individual_method_array, uncertainty_flags, 0.25*intensity_threshold, ...
            %                     [window_size_x, window_size_y], [window_resolution_x, window_resolution_y], display_resampled_images);

            [unc_resampling{trial_index}, snr_metric_resampling{trial_index}, U_resampling, V_resampling] = calculate_uncertainty_resampling(im1_sub, im2_sub, x1, y1, x2, y2, ...
                                d_p*ones(size(x1)), num_particles_remove, num_resampling_trials, individual_method_array, 0.25*intensity_threshold, ...
                                current_pass_settings, display_resampled_images);

            % calculate errors for the resampled uncertainties
            err_U_resampling{trial_index} = U_resampling + results.Ubs(r_trial, c_trial) - errors.Ut_interp(r_trial, c_trial, snapshot_index);
            err_V_resampling{trial_index} = V_resampling + results.Vbs(r_trial, c_trial) - errors.Vt_interp(r_trial, c_trial, snapshot_index);
            
            % ================================================
            %% store error and uncertainty
            % ================================================
            % from original processing
            unc_trials{trial_index} = unc_current;
            snr_metric_trials{trial_index} = snr_metric_current;
            err_U(trial_index) = err_U_current;
            err_V(trial_index) = err_V_current;

            % from windowed processing
            unc_sub_trials{trial_index} = unc_sub_current;
            snr_metric_sub_trials{trial_index} = snr_metric_sub_current;
            err_U_sub_trials(trial_index) = err_U_sub_current;
            err_V_sub_trials(trial_index) = err_V_sub_current;
             
        end

        % ================================================
        %% save calculations to file
        % ================================================
        fprintf('saving results to file\n');
        filename = fullfile(current_write_directory, 'monte_carlo_results.mat');
        save(filename, 'err_U', 'err_V', 'unc_trials', 'snr_metric_trials', ...
            'unc_resampling', 'snr_metric_resampling', 'err_U_resampling', 'err_V_resampling', ...
            'unc_sub_trials', 'snr_metric_sub_trials', 'err_U_sub_trials', 'err_V_sub_trials');
    end
end
% stop timer
toc

% shut down parallel pool
shut_down_parpool();
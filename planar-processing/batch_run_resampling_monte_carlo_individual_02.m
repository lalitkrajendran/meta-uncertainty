% This script performs resampling calculations for a range of particle removal
% percentage and plots the change in uncertainties for different methods

clear
close all
clc

restoredefaultpath;
addpath(genpath('/scratch/shannon/c/aether/Projects/BOS/general-codes/matlab-codes/'));
addpath(genpath('/scratch/shannon/c/aether/Projects/BOS/general-codes/piv-image-generation/'));
addpath('../prana/')
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
dataset_name_array = {'PivChal03B'}; % 'PivChal05B'; 'stagnation_flow'; 'Vortex_Ring'; 'Jetdata'};
num_datasets = numel(dataset_name_array);
% array of uncertainty methods
individual_method_array = {'IM'; 'MC'; 'CS'};
num_individual_methods = numel(individual_method_array);
% array of snr metrics
snr_metric_array = {'PPR'; 'MI'};
num_snr_methods = numel(snr_metric_array);
% directory where results of this analysis are to be saved
top_write_directory = fullfile('/scratch/shannon/a/aether/Projects/PIV_PTV_Uncertainty_Quantification/uncertainty-combination/analysis/results/planar-dataset', 'monte-carlo', 'individual-datasets');
mkdir_c(top_write_directory);
% number of trials
num_trials = 1; %1e3;

% ============================
%% resampling settings
% ============================
% resampling method names
resampling_method_name_array = {'remove-paired'; 'remove-random'; 'add-random'};
% resampling_method_name_array = {'add-random'};
num_resampling_methods = numel(resampling_method_name_array);
% number of particles to remove
percentage_particles_remove_array = 0.05:0.05:0.25; %0.025:0.025:0.25; 
num_ppr = numel(percentage_particles_remove_array);
% number of resampling trials
num_resampling_trials = 1e3; %1e2;
% display_resampled_images? (true/false)
display_resampled_images = 0;
% % case name
% resampling_case_name = ['ppr' num2str(percentage_particles_remove, '%.2f') '_nrs' num2str(num_resampling_trials) '-prewindow-randpart-snr-newproc'];

% ============================
%% particle identification settings
% ============================
% intensity threshold for particle identification
intensity_threshold_percentile = 10;
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
%% plot settings
% ========================
% display figures?
display_figures = 1; %0;
% user screen resolution
user_screen_resolution = 113;
% save_figure? (true/false)
save_figures = 1; %0; %1;
% symbols
symbols = {'o', '^', 'v'};
% colors
colors = lines(3);
% line symbols
line_symbols = {'-'; ':'; '-.'};
% range of displacements to be displayed in the contour plots
displacement_color_min = [0, 0, -0.25, -2, 0];
displacement_color_max = [15, 5, 0.25, 2, 5]; 
% displacement_contour_levels = linspace(displacement_color_min, displacement_color_max, 100);

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

fprintf('running monte-carlo\n');
% set seed for random number generator
% rng(0);
rng();

snapshot_all = cell(num_window_resolution, num_datasets);
X_all = cell(num_window_resolution, num_datasets);
Y_all = cell(num_window_resolution, num_datasets);
unc_mean_all = cell(num_window_resolution, num_datasets);
d_unc_mean_all = cell(num_window_resolution, num_datasets);
snr_mean_all = cell(num_window_resolution, num_datasets);

% start timer
tic
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

        % directory to save results for this case
        % current_write_directory = fullfile(top_write_directory, dataset_name, ['WS' num2str(window_resolution_index)], 'resampling-study');
        current_write_directory = fullfile(top_write_directory, dataset_name, ['WS' num2str(window_resolution_index)], 'weights-rms-change-study');        
        mkdir_c(current_write_directory);

        % load errors
        errors = errors_all{window_resolution_index, dataset_index};

        % calculate random number for data set        
        dataset_index = dataset_index;

        % calculate number of snapshots available for the current dataset and
        % window size
        num_snapshots = numel(results_all{window_resolution_index, dataset_index});
        if num_snapshots > size(errors.err_U, 3)
            num_snapshots = size(errors.err_U, 3);
        end

        % ====================================
        % extract job file properties
        % ====================================
        % extract job file
        jobfile = jobfile_all{window_resolution_index, dataset_index}.Data;
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
            % extract region with non-nan values
            mask = ~isnan(errors.err_U(:, :, 1));
            [cmin, cmax, rmin, rmax] = extract_image_mask_extents(mask(:, :, 1));
            
            % % extract co-ordinate arrays for the processing grid
            % xvec = errors.X(1, :);
            % yvec = errors.Y(:, 1);
            
            % find vector locations in the processing results for which
            % true solution is available
            r = rmin:rmax;
            c = cmin:cmax;
            % c = find(xvec > errors.X(1, 1) & xvec < errors.X(1, end));
            % r = find(yvec > errors.Y(1, 1) & yvec < errors.Y(end, 1));
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
        window_resolution_index_array = nans(1, num_trials);
        
        % uncertainties
        unc_sub_trials = cell(1, num_trials);
        unc_resampling_trials = cell(num_trials, num_ppr, 3);

        % snr metric
        snr_sub_trials = cell(1, num_trials);
        snr_resampling_trials = cell(num_trials, num_ppr, 3);                        
        
        % locations
        snapshot_trials = nans(1, num_trials);
        r_trials = nans(1, num_trials);
        c_trials = nans(1, num_trials);
        X_trials = nans(1, num_trials);
        Y_trials = nans(1, num_trials);
        X_all{window_resolution_index, dataset_index} = nans(1, num_trials);
        Y_all{window_resolution_index, dataset_index} = nans(1, num_trials);
        snapshot_all{window_resolution_index, dataset_index} = nans(1, num_trials);        
        unc_trials = cell(1, num_trials);

        % ================================================
        %% loop through trials
        % ================================================
        for trial_index = 1:num_trials
            fprintf('trial_index: %d\n', trial_index);
            % record dataset index and window resolution index
            window_resolution_index_array(trial_index) = window_resolution_index;
            dataset_index_array(trial_index) = dataset_index;
            
            if window_resolution_index == 1
                % calculate random number for the snapshot
                snapshot_index = randi(num_snapshots);
                % selectr grid point
                r_trial = randsample(r, 1);
                c_trial = randsample(c, 1);
            else
                % for the second window resolution, ensure that the same (or closest) grid location has been selected as for the 1st
                snapshot_index = snapshot_all{1, dataset_index};
                [~, r_trial] = min(abs(results_all{window_resolution_index, dataset_index}{snapshot_index}.Y(:, 1) - Y_all{1, dataset_index}(trial_index)));
                [~, c_trial] = min(abs(results_all{window_resolution_index, dataset_index}{snapshot_index}.X(1, :) - X_all{1, dataset_index}(trial_index))); 
            end

            % load dataset
            results = results_all{window_resolution_index, dataset_index}{snapshot_index};

            % calculate size of the results array
            [num_rows, num_cols] = size(results.X);

            % ================================================
            %% extract processing results for this grid point
            % ================================================
            % extract co-ordinates
            X = results.X(r_trial, c_trial);
            Y = results.Y(r_trial, c_trial);
            
            % store snapshot
            snapshot_all{window_resolution_index, dataset_index} = snapshot_index;

            % store co-ordinates
            X_all{window_resolution_index, dataset_index}(trial_index) = X;
            Y_all{window_resolution_index, dataset_index}(trial_index) = Y;

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
                r_t = find(errors.Y(:, 1) == Y);
                c_t = find(errors.X(1, :) == X);
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
            % im1 = imread(fullfile(files_im1{window_resolution_index, dataset_index}(snapshot_index).folder, files_im1{window_resolution_index, dataset_index}(snapshot_index).name));
            % im1d = double(flipud(im1));
            saved_img = load(fullfile(files_im1{window_resolution_index, dataset_index}(snapshot_index).folder, files_im1{window_resolution_index, dataset_index}(snapshot_index).name));
            im1d = saved_img.im1d;
                
            % im2 = imread(fullfile(files_im2{window_resolution_index, dataset_index}(snapshot_index).folder, files_im2{window_resolution_index, dataset_index}(snapshot_index).name));
            % im2d = double(flipud(im2));
            saved_img = load(fullfile(files_im2{window_resolution_index, dataset_index}(snapshot_index).folder, files_im2{window_resolution_index, dataset_index}(snapshot_index).name));        
            im2d = saved_img.im2d;
                
            % extract image height and width
            if trial_index == 1
                [image_height, image_width] = size(im1d);
            end

            % extract interrogation window
            im1_sub = extract_interrogation_window(im1d, X, Y, [window_size_x, window_size_y]);
            im2_sub = extract_interrogation_window(im2d, X, Y, [window_size_x, window_size_y]);

            % calculate uncertainties
            [unc_sub_trials{trial_index}, snr_sub_trials{trial_index}] = calculate_planar_uncertainties(im1_sub, im2_sub, current_pass_settings, individual_method_array);

            % calculate average particle diameter
            d_p = results.uncertainty2D.Autod(r_trial, c_trial)/sqrt(2);

            % ================================================
            %% loop through particle removal percentages
            % ================================================
            for particle_remove_index = 1:num_ppr
                % current percentage to remove
                percentage_particles_remove = percentage_particles_remove_array(particle_remove_index);
                fprintf('particle remove percentage: %.2f\n', percentage_particles_remove);
                
                % ==========================
                % loop through resampling methods
                % ==========================
                for resampling_method_index = 1:num_resampling_methods
                    fprintf('resampling method index: %d\n', resampling_method_index);
                    resampling_method_name = resampling_method_name_array{resampling_method_index};
                    % calculate uncertainty from resampling
                    [unc, snr] = calculate_uncertainty_resampling_general(im1_sub, im2_sub, d_p, sizeprops, resampling_method_name, percentage_particles_remove, ...
                                                                            num_resampling_trials, individual_method_array, intensity_threshold_percentile, ...
                                                                            current_pass_settings, display_resampled_images);
                    
                    % save_figure_to_png_svg_fig(current_write_directory, ['perturbed-image-' num2str(resampling_method_index)], [1, 0, 0]);
                    % save_figure_to_png_svg_fig('./', 'resampled-image-zoomed-out', [1, 0, 0]);
                    % return;
                    % copy results
                    unc_resampling_trials{trial_index, particle_remove_index, resampling_method_index} = unc;
                    snr_resampling_trials{trial_index, particle_remove_index, resampling_method_index} = snr;
                end    

            end

            % ================================================
            %% store uncertainties and locations
            % ================================================
            snapshot_trials(trial_index) = snapshot_index;
            r_trials(trial_index) = r_trial;
            c_trials(trial_index) = c_trial;            
            X_trials(trial_index) = X;
            Y_trials(trial_index) = Y;            
            unc_trials{trial_index} = unc_current;            
        end

        % ================================================
        % save results
        % ================================================
        
        % file name
        filename = ['resampling-statistics-nt=' num2str(num_trials, '%d') 'nr='  num2str(num_resampling_trials, '%d') '.mat'];
        % save results to file
        save(fullfile(current_write_directory, filename), 'unc_sub_trials', 'snr_sub_trials', 'unc_trials', ...
                        'unc_resampling_trials', 'snr_resampling_trials', ...
                        'r_trials', 'c_trials', 'snapshot_trials', 'X_trials', 'Y_trials', ...
                        'percentage_particles_remove_array', 'resampling_method_name_array');   
        break;                      
    end
end

% shut down parallel pool
shut_down_parpool();

% stop timer
toc

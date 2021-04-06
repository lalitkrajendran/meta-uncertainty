% This script performs resampling calculations for a range of particle removal
% percentage and plots the change in uncertainties for different methods

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
num_individual_methods = numel(individual_method_array);
% array of snr metrics
snr_metric_array = {'PPR'; 'MI'};
num_snr_methods = numel(snr_metric_array);
% directory where results of this analysis are to be saved
top_write_directory = fullfile('/scratch/shannon/a/aether/Projects/PIV_PTV_Uncertainty_Quantification/uncertainty-combination/analysis/results/planar-dataset', 'monte-carlo', 'individual-datasets');
mkdir_c(top_write_directory);
% number of trials
num_trials = 1e1;
% displacement component names
component_names = {'x'; 'y'};
num_components = numel(component_names);

% ============================
%% resampling settings
% ============================
% uncertainty deviation percentage
uncertainty_deviation_percentage = 0.1;
% number of particles to remove
% num_particles_remove_array = [5, 10, 15];
percentage_particles_remove = 0.1;
% percentage_particles_remove_array = [0.1]; %0:0.025:0.25;
percentage_particles_remove_array = 0:0.025:0.10; %0.25;
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
for window_size_index = 1:num_window_resolution
    % loop through datasets
    for dataset_index = 1:num_datasets
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
    for dataset_index = 1:num_datasets
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
        ['trials=' num2str(num_trials, '%d')]);
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
        
        % locations
        r_trials = nans(1, num_trials);
        c_trials = nans(1, num_trials);
        X_trials = nans(1, num_trials);
        Y_trials = nans(1, num_trials);
        X_all{window_size_index, dataset_index} = nans(1, num_trials);
        Y_all{window_size_index, dataset_index} = nans(1, num_trials);
        
        % ================================================
        % load results
        % ================================================
        save_directory = fullfile(top_write_directory, dataset_name, ['WS' num2str(window_size_index)], ['trials=' num2str(num_trials, '%d')]);
        % resampling_statistics = load(fullfile(save_directory, 'resampling-statistics.mat'));
        resampling_statistics = load(fullfile(save_directory, ['resampling-statistics-' num2str(num_resampling_trials, '%d') '.mat']));

        ppr_thresh = nans(num_trials, num_individual_methods, num_components);
        min_ppr_thresh = nans(num_trials, 1);
        % calculate minimum allowable threshold for specified uncertainty deviation
        for trial_index = 1:num_trials
            % calculate random number for the snapshot
            snapshot_index = 1; %randi(num_snapshots);

            % load dataset
            results = results_all{window_size_index, dataset_index}{snapshot_index};

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

            fprintf('trial_index: %d\n', trial_index);
            % record dataset index and window resolution index
            window_size_index_array(trial_index) = window_size_index;
            dataset_index_array(trial_index) = dataset_index;

            % calculate random number for the snapshot
            snapshot_index = 1; %randi(num_snapshots);

            % load dataset
            results = results_all{window_size_index, dataset_index}{snapshot_index};

            % calculate size of the results array
            [num_rows, num_cols] = size(results.X);
            
            % selectr grid point
            if window_size_index == 1
                r_trial = randsample(r, 1);
                c_trial = randsample(c, 1);
            else
                [~, r_trial] = min(abs(results.Y(:, 1) - Y_all{1, dataset_index}(trial_index)));
                [~, c_trial] = min(abs(results.X(1, :) - X_all{1, dataset_index}(trial_index))); 
            end

            % ================================================
            %% extract processing results for this grid point
            % ================================================
            % extract co-ordinates
            X = results.X(r_trial, c_trial);
            Y = results.Y(r_trial, c_trial);
            
            % store co-ordinates
            X_all{window_size_index, dataset_index}(trial_index) = X;
            Y_all{window_size_index, dataset_index}(trial_index) = Y;

            % ==========================
            %% calculate particle removal percentage threshold
            % ==========================
            for individual_method_index = 1:num_individual_methods
                individual_method_name = lower(individual_method_array{individual_method_index});
                for component_index = 1:num_components
                    % extract uncertainty for this component
                    unc = resampling_statistics.unc.([individual_method_name component_names{component_index}])(trial_index, :);
                    % only retain unique values
                    [u, i] = unique(unc, 'stable');
                    % calculate threshold for specified uncertainty deviation
                    % ppr = interp1(abs(u - u(1))/u(1), percentage_particles_remove_array(i), uncertainty_deviation_percentage, 'pchip', 'extrap');
                    ppr = interp1(abs(u - u(1))/u(1), percentage_particles_remove_array(i), uncertainty_deviation_percentage, 'pchip'); %, 'extrap');
                    ppr_thresh(trial_index, individual_method_index, component_index) = ppr;
                    if ppr < 0
                        fprintf('trial: %d, method: %d, component: %d, ppr: %.2f\n', trial_index, individual_method_index, component_index, ppr*100);
                    end
                end                
            end

            % calculate minimum across all schemes and components
            for component_index = 1:num_components
                ppr = squeeze(ppr_thresh(trial_index, :, component_index));
                [minval, minloc] = min(ppr);
                min_ppr_thresh(trial_index, component_index) = minval;
                min_method(trial_index, component_index) = minloc;
            end
        end
        
        if display_figures
            % ================================================
            %% plot change in uncertainties (contour)
            % ================================================
            figure
            for individual_method_index = 1:num_individual_methods
                individual_method_name = lower(individual_method_array{individual_method_index});
                subplot(num_individual_methods, 1, individual_method_index)
                displacement_contour_levels = linspace(0, displacement_color_max(dataset_index), 100);
                % plot displacements
                contourf(results.X, results.Y, abs(results.U), displacement_contour_levels, 'edgecolor', 'none')
                % colormap(gray)
                % colormap(summer)
                cmap = cbrewer('seq', 'Greys', 100);        
                colormap(cmap);
                caxis([0 displacement_color_max(dataset_index)])
                colorbar
                set_axes(gca)
                axis off    
                hold on
        
                % plot changes in uncertainty
                % d_unc = diff(unc_mean_trials.([individual_method_name 'x']), 1, 2);                
                ppr = squeeze(ppr_thresh(:, individual_method_index, 1)) * 100;                
                plot_colored_circles_seq(gcf, X_all{window_size_index, dataset_index}', Y_all{window_size_index, dataset_index}', ppr, ...
                                            0, 10, min(image_height, image_width)/20);                                    
                title(upper(individual_method_name))            
            end

            % plot circles to highlight the method with minimum threshold
            for trial_index = 1:num_trials
                subplot_index = min_method(trial_index, 1);
                x = X_all{window_size_index, dataset_index}(trial_index);
                y = Y_all{window_size_index, dataset_index}(trial_index);
                subplot(num_individual_methods, 1, subplot_index)                
                plot(x, y, 'o', 'markersize', min(image_height, image_width)/20 * 0.5, 'color', 'g')
            end
            
            sgtitle('Range: 0 - 10 %', 'fontweight', 'bold')
            set(gcf, 'resize', 'off')
            set(gcf, 'units', 'inches', 'position', [381    50   784   684]/user_screen_resolution)
            drawnow();

            if save_figures            
                save_figure_to_png_svg_fig(save_directory, ['particle-removal-threshold-contour-' num2str(num_resampling_trials, '%d')], [1, 0, 0]);
            end            
        end
    end
end

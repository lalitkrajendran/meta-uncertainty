clear
close all
clc

restoredefaultpath;
addpath(genpath('/scratch/shannon/c/aether/Projects/BOS/general-codes/matlab-codes/'));

%% processing settings

% window resolution
window_size_array = [32, 64];
num_window_resolution = numel(window_size_array);

% directory containing images
top_image_directory = fullfile('/scratch/shannon/a/aether/Projects/PIV_PTV_Uncertainty_Quantification/Planar_uncertainty_work/Images/experiment/');

% array of data set names
dataset_name_array = {'PivChal03B'; 'PivChal05B'; 'stagnation_flow'; 'Vortex_Ring'; 'Jetdata'};
num_datasets = numel(dataset_name_array);

% directory containing sample job files
sample_job_file_directory = '/scratch/shannon/a/aether/Projects/PIV_PTV_Uncertainty_Quantification/uncertainty-combination/analysis/data/sample-job-files/';

% directory where results of this analysis are to be saved
top_write_directory = fullfile('/scratch/shannon/a/aether/Projects/PIV_PTV_Uncertainty_Quantification/uncertainty-combination/analysis/results/', 'experiment-new');

% pass number
pass_number = 4;

% directory containing true solution
true_solution_directory = fullfile('/scratch/shannon/a/aether/Projects/PIV_PTV_Uncertainty_Quantification/Planar_uncertainty_work/Results/07_30_2018/True_solution');

%% plot settings

% save figures? (True/False)
save_figures = true;

% minimum uncertainty for a measurement to be considered valid
min_error_threshold = 1e-4;

% multiplication factor for quiver plot
displacement_scale = 10;

% range of displacements to be displayed in the contour plots
displacement_color_min = 0;
displacement_color_max = 15;
displacement_contour_levels = linspace(displacement_color_min, displacement_color_max, 100);

% range of errors to be displayed in the contour plots
err_color_min = 0;
err_color_max = 0.25;
err_contour_levels = linspace(err_color_min, err_color_max, 100);

%% load and display results

for window_size_index = 1:numel(window_size_array)
    fprintf('WS: %d\n', window_size_index);
        
    for dataset_index = 4 %1:3 %1:num_datasets
        % extract data set name
        dataset_name = dataset_name_array{dataset_index};
        fprintf('dataset name: %s\n', dataset_name);
        
        % image directory for current data set
        current_image_directory = fullfile(top_image_directory, dataset_name);
        
        % results directory for current data set
        current_results_directory = fullfile(top_write_directory, dataset_name, ['WS' num2str(window_size_index)]);

        % directory containing vectors
        vectors_directory = fullfile(current_results_directory, 'vectors');

        % get list of mat files containing results for the desired pass
        [files, num_files] = get_directory_listing(vectors_directory, ['*pass' num2str(pass_number) '*.mat']);
        
        % load true solution
        true_solution = load(fullfile(true_solution_directory, dataset_name, 'true_solution.mat'));

        %% initialize error matrices
        
        % load sample result
        sample_results = load(fullfile(files(1).folder, files(1).name));
        % calculate size of arrays
        [num_rows, num_cols] = size(sample_results.X);
        
        % initialize error arrays
        err_U = nans(num_rows, num_cols, num_files);
        err_V = nans(num_rows, num_cols, num_files);
        
        %% prepare true solution 

        % vortex ring
        if dataset_index == 4
            % convert true solution to the image space
            [true_solution, c, r] = calculate_true_solution_vortex_ring(true_solution, sample_results.X, sample_results.Y);
            % extract size of the points where the error is to be
            % calculated
            num_rows_ind = numel(r);
            num_cols_ind = numel(c);
            % initialize co-ordinate arrays
            X_err = sample_results.X(r, c);
            Y_err = sample_results.Y(r, c);
            
            % initialize velocity arrays
            U = nans(num_rows_ind, num_cols_ind, num_files);
            V = nans(num_rows_ind, num_cols_ind, num_files);
            Ut = nans(num_rows_ind, num_cols_ind, num_files);
            Vt = nans(num_rows_ind, num_cols_ind, num_files);
            
            % initialize error arrays
            err_U = nans(num_rows_ind, num_cols_ind, num_files);
            err_V = nans(num_rows_ind, num_cols_ind, num_files);
        end

        %% loop through files and calculate errors
        for file_index = 1:50 %num_files
            % load result
            results = load(fullfile(files(file_index).folder, files(file_index).name));
            
            % interpolate processed displacements
            U(:, :, file_index) = results.U(r, c);
            V(:, :, file_index) = results.V(r, c);
            % interpolate true solution
            Ut(:, :, file_index) = interp2(true_solution.Xt_scaled, true_solution.Yt_scaled, true_solution.Ut_scaled(:, :, file_index), ...
                X_err, Y_err, 'linear', 0);
            Vt(:, :, file_index) = interp2(true_solution.Xt_scaled, true_solution.Yt_scaled, true_solution.Vt_scaled(:, :, file_index), ...
                X_err, Y_err, 'linear', 0);
        end
        
        %% contour plots of error
        file_index = 25;
        
        figure
        imagesc(U(:, :, file_index) - Ut(:, :, file_index))
        caxis([-0.5 0.5])
        title('0')
        
        figure
        imagesc(U(:, :, file_index + 3) - Ut(:, :, file_index))
        caxis([-0.5 0.5])
        title('+3')

        figure
        imagesc(U(:, :, file_index - 3) - Ut(:, :, file_index))
        caxis([-0.5 0.5])
        title('-3')
        
        return;
        %% make quiver plots
        
        figure
        quiver(X_err, Y_err, U(:, :, file_index) * displacement_scale, V(:, :, file_index) * displacement_scale, 'AutoScale', 'off')
        hold on
        quiver(X_err, Y_err, Ut(:, :, file_index) * displacement_scale, Vt(:, :, file_index) * displacement_scale, 'AutoScale', 'off')
        
        %% make line plots
        file_index = 25;
        
        r_plot = find(Y_err(:, 1) == 616);
        c_plot = find(X_err(1, :) == 496);

        figure
        plot(X_err(1,:), Ut(r_plot, :, file_index), 'o')
        hold on
        plot(X_err(1,:), U(r_plot, :, file_index-5), '*')
        plot(X_err(1,:), U(r_plot, :, file_index), '*')
        plot(X_err(1,:), U(r_plot, :, file_index+5), '*')
        
        legend('True', '-5', '0', '+5')
        return;
    end
end
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
        
    for dataset_index = 1:num_datasets
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
        
        % load sample result
        sample_results = load(fullfile(files(1).folder, files(1).name));
        % calculate size of arrays
        [num_rows, num_cols] = size(sample_results.X);

        %% load true solution
        
        % load true solution
        true_solution = load(fullfile(true_solution_directory, dataset_name, 'true_solution.mat'));
        
        % set the number of files to be the lowest between the true
        % solution and the processing
        if dataset_index == 4
            num_files = 50;
        elseif dataset_index == 5
            num_files = 496;
        end        
        
        %% prepare true solution 
        
        if dataset_index < 4
            num_rows_p = num_rows;
            num_cols_p = num_cols;
            X = sample_results.X;
            Y = sample_results.Y;
        % vortex ring
        elseif dataset_index == 4
            % convert true solution to the image space
            [true_solution, c, r] = calculate_true_solution_vortex_ring(true_solution, sample_results.X, sample_results.Y);
            % extract size of the points where the error is to be
            % calculated
            num_rows_p = numel(r);
            num_cols_p = numel(c);

            % initialize co-ordinate arrays
            X = sample_results.X(r, c);
            Y = sample_results.Y(r, c);
            
        elseif dataset_index == 5
            % index limits for processing
            % X 320 to 392 (72); Y 205 to 325 (120)
            cmin_p = 5;
            cmax_p = 23; %303+18-1=320  303+90-1=392  X direction
            rmin_p = 5;
            rmax_p = 35; %188+18-1=205  188+138-1=325 Y direction
            
            % number of indices for processing
            num_rows_p = rmax_p - rmin_p + 1;
            num_cols_p = cmax_p - cmin_p + 1;
            
            % indices to be used
            rp = rmin_p:rmax_p;
            cp = cmin_p:cmax_p;
            
            % indices to be used in the true solution
            rt = 8:2:69;
            ct = 7:2:43;
            
            % extract co-ordinate arrays
            X = sample_results.X(rp, cp) + 302;
            Y = sample_results.Y(rp, cp) + 187;
        end
        
        %% initialize arrays
        
        % initialize velocity arrays
        U = nans(num_rows_p, num_cols_p, num_files);
        V = nans(num_rows_p, num_cols_p, num_files);
        Ut = nans(num_rows_p, num_cols_p, num_files);
        Vt = nans(num_rows_p, num_cols_p, num_files);
        
        % initialize error arrays
        err_U = nans(num_rows_p, num_cols_p, num_files);
        err_V = nans(num_rows_p, num_cols_p, num_files);

        %% loop through files and calculate errors
        
        for file_index = 1:num_files
            % load result
            if dataset_index < 5
                results = load(fullfile(files(file_index).folder, files(file_index).name));
            else
                results = load(fullfile(files(file_index+1).folder, files(file_index+1).name));
            end
                        
            % calculate error
            if dataset_index == 1
                err_U(:, :, file_index) = results.U - rot90(true_solution.u(:, :, file_index));
                err_V(:, :, file_index) = results.V - rot90(true_solution.v(:, :, file_index));
            elseif dataset_index == 2
                err_U(:, :, file_index) = results.U - true_solution.Ut(:, :, file_index);
                err_V(:, :, file_index) = results.V - true_solution.Vt(:, :, file_index);                
            elseif dataset_index == 3
                err_U(:, :, file_index) = results.U - true_solution.Ut;
                err_V(:, :, file_index) = results.V - true_solution.Vt;                
            elseif dataset_index == 4
                % interpolate processed displacements
                U(:, :, file_index) = results.U(r, c);
                V(:, :, file_index) = results.V(r, c);
                % interpolate true solution
                Ut(:, :, file_index) = interp2(true_solution.Xt_scaled, true_solution.Yt_scaled, true_solution.Ut_scaled(:, :, file_index), ...
                    X, Y, 'linear', 0);
                Vt(:, :, file_index) = interp2(true_solution.Xt_scaled, true_solution.Yt_scaled, true_solution.Vt_scaled(:, :, file_index), ...
                    X, Y, 'linear', 0);
                % calculate error
                err_U(:, :, file_index) = U(:, :, file_index) - Ut(:, :, file_index);
                err_V(:, :, file_index) = V(:, :, file_index) - Vt(:, :, file_index);
            elseif dataset_index == 5
                % calculate error
                err_U(:, :, file_index) = results.U(rp, cp) - flipud(true_solution.Ut(rt, ct, file_index));
                err_V(:, :, file_index) = results.V(rp, cp) - flipud(true_solution.Vt(rt, ct, file_index));                
            end
            
        end
        
        %% calculate error statistics
        
        % U
        [bias_U, random_U, total_U, rms_U] = calculate_error_statistics(err_U, 3);
        % V
        [bias_V, random_V, total_V, rms_V] = calculate_error_statistics(err_V, 3);
        
        % save errors to file
        save(fullfile(current_results_directory, 'errors.mat'), 'X', 'Y', 'err_U', 'err_V', ...
                                                'bias_U', 'random_U', 'total_U', ...
                                                'bias_V', 'random_V', 'total_V');

        %% display error statistics
        
        % directory to save figures
        figure_write_directory = fullfile(current_results_directory, 'figures-axis-off');
        mkdir_c(figure_write_directory);
        
        if dataset_index < 4
            U_plot = sample_results.U;
            V_plot = sample_results.V;
        elseif dataset_index == 4
            U_plot = U(:, :, 1);
            V_plot = V(:, :, 1);            
        elseif dataset_index == 5
            U_plot = sample_results.U(rp, cp);
            V_plot = sample_results.V(rp, cp);            
        end
        figure        
        % displacement
        subplot(2, 2, 1)
        contourf(X, Y, sqrt(U_plot.^2 + V_plot.^2), displacement_contour_levels, 'edgecolor', 'none');
        colormap(flipud(gray))
        annotate_image(gcf, gca);
        h2 = colorbar;
        caxis([displacement_color_min displacement_color_max])
        axis off
        title(h2, '\Delta x (pix.)')
        title('Displacement')
        
        % bias error
        subplot(2, 2, 2)
        contourf(X, Y, sqrt(bias_U.^2 + bias_V.^2), err_contour_levels, 'edgecolor', 'none');
        colormap(flipud(gray))
        annotate_image(gcf, gca);
        h2 = colorbar;
        caxis([err_color_min err_color_max])
        axis off
        title(h2, '(pix.)')
        title({'Bias'; 'Error'})

        % random error
        subplot(2, 2, 3)
        contourf(X, Y, sqrt(random_U.^2 + random_V.^2), err_contour_levels, 'edgecolor', 'none');
        colormap(flipud(gray))
        annotate_image(gcf, gca);
        h2 = colorbar;
        caxis([err_color_min err_color_max])
        axis off
        title(h2, '(pix.)')
        title({'Random'; 'Error'})

        % total error
        subplot(2, 2, 4)
        contourf(X, Y, sqrt(total_U.^2 + total_V.^2), err_contour_levels, 'edgecolor', 'none');
        colormap(flipud(gray))
        annotate_image(gcf, gca);
        h2 = colorbar;
        caxis([err_color_min err_color_max])
        axis off
        title(h2, '(pix.)')
        title({'Total'; 'Error'})
        
        % adjust figure size for each dataset
        if dataset_index == 1
            set(gcf, 'Position', [56   285   1204   588]);
        elseif dataset_index == 2
            set(gcf, 'Position', [56   335   1147   538]);
        elseif dataset_index == 3
            set(gcf, 'Position', [167   186   895   676]);
        elseif dataset_index == 4
            set(gcf, 'Position', [173   263   953   636]);            
        elseif dataset_index == 5
            set(gcf, 'Position', [150   240   835   640]);                    
        end

        % save figures if needed. else pause to view the figures
        if save_figures
            save_figure_to_png_eps_fig(figure_write_directory, 'error-statistics-contours', [true, false, false]);
        else
            pause(0.1);
        end
        % close figure
        close(gcf);

    end
end
clear
close all
clc

restoredefaultpath;
addpath(genpath('/scratch/shannon/c/aether/Projects/BOS/general-codes/matlab-codes/'));
addpath('../general-codes');
setup_default_settings;

% ====================
%% processing settings
% ====================
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
top_write_directory = fullfile('/scratch/shannon/a/aether/Projects/PIV_PTV_Uncertainty_Quantification/uncertainty-combination/analysis/results/planar-dataset', 'experiment-new');
% pass number
pass_number = 4;
% directory containing true solution
true_solution_directory = fullfile('/scratch/shannon/a/aether/Projects/PIV_PTV_Uncertainty_Quantification/Planar_uncertainty_work/Results/07_30_2018/True_solution');

% ====================
%% plot settings
% ====================
% display figures? (True/False)
display_figures = 1;
% save figures? (True/False)
save_figures = 0;
% user screen resolution
user_screen_resolution = 113;
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

% ====================
%% load and display results
% ====================
for window_size_index = 1:numel(window_size_array)
    fprintf('WS: %d\n', window_size_index);        
    for dataset_index = 5 %1:num_datasets        
        % extract data set name
        dataset_name = dataset_name_array{dataset_index};
        fprintf('dataset name: %s\n', dataset_name);
        
        % image directory for current data set
        current_image_directory = fullfile(top_image_directory, dataset_name);
        
        % results directory for current data set
        current_results_directory = fullfile(top_write_directory, dataset_name, ['WS' num2str(window_size_index)]);

        % directory containing vectors
        vectors_directory = fullfile(current_results_directory, 'vectors-new');

        % get list of mat files containing results for the desired pass
        [files, num_files] = get_directory_listing(vectors_directory, ['*pass' num2str(pass_number) '*.mat']);
        
        % load sample result
        sample_results = load(fullfile(files(1).folder, files(1).name));
        % calculate size of arrays
        [num_rows, num_cols] = size(sample_results.X);

        % ====================
        %% load true solution
        % ====================        
        % load true solution
        true_solution = load(fullfile(true_solution_directory, dataset_name, 'true_solution.mat'));
        
        % set the number of files to be the lowest between the true
        % solution and the processing
        if dataset_index == 4
            num_files = 50;
        elseif dataset_index == 5
            num_files = 496;
        end        
        
        % ====================
        %% prepare true solution 
        % ====================        
        if dataset_index < 4
            %% piv challenge 03B and 05B and stagnation flow
            num_rows_p = num_rows;
            num_cols_p = num_cols;
            X = sample_results.X;
            Y = sample_results.Y;
        elseif dataset_index == 4
            %% vortex ring
            % convert true solution to the image space
            [true_solution, c, r] = calculate_true_solution_vortex_ring(true_solution, sample_results.X, sample_results.Y);
            % extract size of the points where the error is to be
            % calculated
            num_rows_p = numel(r);
            num_cols_p = numel(c);

            % initialize co-ordinate arrays
            % X = sample_results.X(r, c);
            % Y = sample_results.Y(r, c);
            X = sample_results.X;
            Y = sample_results.Y;
            
        elseif dataset_index == 5
            %% jet
            % index limits for processing
            % X 320 to 392 (72); Y 205 to 325 (120)
            % extract co-ordinate arrays
            X = sample_results.X + 302;
            Y = sample_results.Y + 187;
        end
        
        % ====================
        %% loop through files and calculate errors
        % ====================
        
        % initialize error arrays
        Ut_interp = nans(size(sample_results.X, 1), size(sample_results.X, 2), num_files);
        Vt_interp = nans(size(sample_results.X, 1), size(sample_results.X, 2), num_files);
        err_U = nans(size(sample_results.X, 1), size(sample_results.X, 2), num_files);
        err_V = nans(size(sample_results.X, 1), size(sample_results.X, 2), num_files);

        % directory to save figures
        if save_figures
            figure_write_directory = fullfile(current_results_directory, 'figures-new');
            mkdir_c(figure_write_directory);
        end

        for file_index = 1:num_files
            % ====================
            %% load result
            % ====================
                        
            if dataset_index < 5
                results = load(fullfile(files(file_index).folder, files(file_index).name));
            else
                results = load(fullfile(files(file_index+1).folder, files(file_index+1).name));
            end
            
            % ====================
            %% modify true solution
            % ====================
            if dataset_index == 1
                Xt = rot90(true_solution.x);
                Yt = rot90(true_solution.y);
                Ut = rot90(true_solution.u(:, :, file_index));
                Vt = rot90(true_solution.v(:, :, file_index));
            elseif dataset_index == 2
                Xt = true_solution.Xt;
                Yt = true_solution.Yt;
                Ut = true_solution.Ut(:, :, file_index);
                Vt = true_solution.Vt(:, :, file_index);
            elseif dataset_index == 3
                xt = 8 + 16 * (0:79);
                yt = 8 + 16 * (0:63);
                [Xt, Yt] = meshgrid(xt, yt);
                Ut = true_solution.Ut;
                Vt = true_solution.Vt;
            elseif dataset_index == 4
                Xt = true_solution.Xt_scaled;
                Yt = true_solution.Yt_scaled;
                Ut = true_solution.Ut_scaled(:, :, file_index);
                Vt = true_solution.Vt_scaled(:, :, file_index);
            elseif dataset_index == 5
                Xt = true_solution.X;
                Yt = true_solution.Y;
                Ut = flipud(true_solution.Ut(:, :, file_index));
                Vt = -flipud(true_solution.Vt(:, :, file_index));
                % only retain true solution in the index limits for
                % processing
                c = Xt(1, :) >= 320  & Xt(1, :) <= 392;
                r = Yt(:, 1) >= 205  & Yt(:, 1) <= 325;
                Xt = Xt(r, c);
                Yt = Yt(r, c);
                Ut = Ut(r, c);
                Vt = Vt(r, c);
            end
            
            % ====================
            %% interpolate true solution onto the measurement grid
            % ====================            
            Ut_interp(:, :, file_index) = interp2(Xt, Yt, Ut, X, Y, 'linear', NaN);
            Vt_interp(:, :, file_index) = interp2(Xt, Yt, Vt, X, Y, 'linear', NaN);

            % ====================
            %% calculate error
            % ====================            
            err_U(:, :, file_index) = results.U - Ut_interp(:, :, file_index);
            err_V(:, :, file_index) = results.V - Vt_interp(:, :, file_index);            
        end
        
        % ====================
        %% plot instantaneous results
        % ====================
        if display_figures
            plot_instantaneous_error(X, Y, sample_results.U, sample_results.V, Ut_interp(:, :, 1), Vt_interp(:, :, 1), err_U(:, :, 1), err_V(:, :, 1));
            
            set(gcf, 'resize', 'off');
            set(gcf, 'units', 'inches', 'position', [35         342        1233         452]/user_screen_resolution);
            drawnow();
        
            if save_figures
                save_figure_to_png_svg_fig(figure_write_directory, 'displacements-errors-instantaneous', [1, 0, 0]);
            end
        end

        % ====================
        %% calculate error statistics
        % ====================        
        % U
        [bias_U, random_U, total_U, rms_U] = calculate_error_statistics(err_U, 3);
        % V
        [bias_V, random_V, total_V, rms_V] = calculate_error_statistics(err_V, 3);
        
        % ====================
        % save errors to file
        % ====================
        save(fullfile(current_results_directory, 'errors-new.mat'), 'X', 'Y', ...
             'Xt', 'Yt', 'Ut_interp', 'Vt_interp', ...
             'err_U', 'err_V', 'bias_U', 'random_U', 'total_U', ...
             'bias_V', 'random_V', 'total_V');
        
        % ====================
        %% display error statistics
        % ====================
        if display_figures        
            U_plot = sample_results.U;
            V_plot = sample_results.V;
            
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
                pos = [56   285   1204   588];
            elseif dataset_index == 2
                pos = [56   335   1147   538];
            elseif dataset_index == 3
                pos = [167   186   895   676];
            elseif dataset_index == 4
                pos = [173   263   953   636];
            elseif dataset_index == 5
                pos = [150   240   835   640];
            end
            set(gcf, 'units', 'inches', 'Position', pos/user_screen_resolution);

            % save figures if needed. else pause to view the figures
            if save_figures
                % save_figure_to_png_eps_fig(figure_write_directory, 'error-statistics-contours', [true, false, false]);
                save_figure_to_png_svg_fig(figure_write_directory, 'error-statistics-contours', [1, 0, 0]);
            else
                pause(0.1);
            end

            % close figure
            close(gcf);
        end
    end
end
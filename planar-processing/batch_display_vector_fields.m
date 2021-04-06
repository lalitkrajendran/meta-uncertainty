clear
close all
clc

restoredefaultpath;
addpath(genpath('/scratch/shannon/c/aether/Projects/BOS/general-codes/matlab-codes/'));
addpath ../prana/

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

%% plot settings

% save figures? (True/False)
save_figures = true;

% minimum uncertainty for a measurement to be considered valid
min_error_threshold = 1e-4;

% range of displacements to be displayed in the contour plots
displacement_color_min = 0;
displacement_color_max = 15;
displacement_contour_levels = linspace(displacement_color_min, displacement_color_max, 100);

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
        results = load(fullfile(files(1).folder, files(1).name));
        
        %% set uncertainty range for the data set
        if dataset_index == 3 || dataset_index == 5
            % maximum error for a measurement to be considered valid
            max_error_threshold = 0.3;
            
            % range of uncertainties to be displayed (pix.)
            uncertainty_color_min = 0;
            uncertainty_color_max = 0.3;
            uncertainty_contour_levels = linspace(uncertainty_color_min, uncertainty_color_max, 100);
            max_bin = 0.3;
        else
            % maximum error for a measurement to be considered valid
            max_error_threshold = 0.1;
            
            % range of uncertainties to be displayed (pix.)
            uncertainty_color_min = 0;
            uncertainty_color_max = 0.1;
            uncertainty_contour_levels = linspace(uncertainty_color_min, uncertainty_color_max, 100);
            max_bin = 0.08;
        end
        %% display result
        
        % directory to save figures
        figure_write_directory = fullfile(current_results_directory, 'figures');
        mkdir_c(figure_write_directory);
        
        figure
        % displacement
        subplot(2, 2, 1)
        contourf(results.X, results.Y, sqrt(results.U.^2 + results.V.^2), displacement_contour_levels, 'edgecolor', 'none');
        colormap(flipud(gray))
        annotate_image(gcf, gca);
        h2 = colorbar;
        caxis([displacement_color_min displacement_color_max])
        title(h2, '\Delta x (pix.)')
        title('Displacement')
        
        % IM
        subplot(2, 2, 2)
        contourf(results.X, results.Y, sqrt(results.uncertainty2D.Uimx.^2 + results.uncertainty2D.Uimy.^2), uncertainty_contour_levels, 'edgecolor', 'none');
        colormap(flipud(gray))
        annotate_image(gcf, gca);
        h2 = colorbar;
        caxis([uncertainty_color_min uncertainty_color_max])
        title(h2, '\sigma_{\Delta x} (pix.)')
        title({'Image'; 'Matching'})

        % MC
        subplot(2, 2, 3)
        contourf(results.X, results.Y, sqrt(results.uncertainty2D.MCx.^2 + results.uncertainty2D.MCy.^2), uncertainty_contour_levels, 'edgecolor', 'none');
        colormap(flipud(gray))
        annotate_image(gcf, gca);
        h2 = colorbar;
        caxis([uncertainty_color_min uncertainty_color_max])
        title(h2, '\sigma_{\Delta x} (pix.)')
        title({'Moment of'; 'Correlation'})

        % CS
        subplot(2, 2, 4)
        contourf(results.X, results.Y, sqrt(real(results.uncertainty2D.Ucsx).^2 + real(results.uncertainty2D.Ucsy).^2), uncertainty_contour_levels, 'edgecolor', 'none');
        colormap(flipud(gray))
        annotate_image(gcf, gca);
        h2 = colorbar;
        caxis([uncertainty_color_min uncertainty_color_max])
        title(h2, '\sigma_{\Delta x} (pix.)')
        title({'Correlation';'Statistics'})
        
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
            save_figure_to_png_eps_fig(figure_write_directory, 'displacements-uncertainties-contours', [true, false, false]);
        else
            pause(0.1);
        end
        % close figure
        close(gcf);
    end
end
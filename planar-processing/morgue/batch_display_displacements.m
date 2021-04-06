clear
close all
clc


restoredefaultpath;
addpath(genpath('/scratch/shannon/c/aether/Projects/BOS/general-codes/matlab-codes/'));
setup_default_settings;
% addpath ../prana/
addpath ../prana-uncertainty-average-dc-new-im-2-save-deform-cs/

%% processing settings

% directory containing images
top_image_directory = fullfile('/scratch/shannon/a/aether/Projects/PIV_PTV_Uncertainty_Quantification/Planar_uncertainty_work/Images/experiment/');

% array of data set names
dataset_name_array = {'PivChal03B'; 'PivChal05B'; 'stagnation_flow'; 'Vortex_Ring'; 'Jetdata'};
num_datasets = numel(dataset_name_array);

% directory containing sample job files
sample_job_file_directory = '/scratch/shannon/a/aether/Projects/PIV_PTV_Uncertainty_Quantification/uncertainty-combination/analysis/data/sample-job-files/';

% directory where results of this analysis are to be saved
top_write_directory = fullfile('/scratch/shannon/a/aether/Projects/PIV_PTV_Uncertainty_Quantification/uncertainty-combination/analysis/results/', 'experiment-new');
if ~exist(top_write_directory, 'dir')
    mkdir(top_write_directory);
end

% pass number
pass_number = 4;

%% plot settings

% save figures? (True/False)
save_figures = false;

% minimum uncertainty for a measurement to be considered valid
min_error_threshold = 1e-4;

% range of displacements to be displayed in the contour plots
displacement_color_min = 0;
displacement_color_max = 15;
displacement_contour_levels = linspace(displacement_color_min, displacement_color_max, 100);

% number of bins for the coarse histogram
num_bins = 8;
%% run processing

for window_size_index = 1:2
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
        
        % load pass results
        results = load_directory_data(vectors_directory, ['*pass' num2str(pass_number) '*.mat']);

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
        %%
        
        % display results
        % =================================
        % Prana, IM and MC contours
        % =================================

        figure
        subplot(2,2,1)
        contourf(results.X, results.Y, sqrt(results.U.^2 + results.V.^2), displacement_contour_levels, 'edgecolor', 'none');
        colormap(flipud(gray))
        annotate_image(gcf, gca);
        h2 = colorbar;
        caxis([displacement_color_min displacement_color_max])
        title(h2, '\Delta x (pix.)')
        title('Displacement')

        subplot(2,2,2)
        contourf(results.X, results.Y, sqrt(results.uncertainty2D.Uimx.^2 + results.uncertainty2D.Uimy.^2), uncertainty_contour_levels, 'edgecolor', 'none');
        colormap(flipud(gray))
        annotate_image(gcf, gca);
        h2 = colorbar;
        caxis([uncertainty_color_min uncertainty_color_max])
        title(h2, '\sigma_{\Delta x} (pix.)')
        title({'Image'; 'Matching'})

        subplot(2,2,3)
        contourf(results.X, results.Y, sqrt(results.uncertainty2D.MCx.^2 + results.uncertainty2D.MCy.^2), uncertainty_contour_levels, 'edgecolor', 'none');
        colormap(flipud(gray))
        annotate_image(gcf, gca);
        h2 = colorbar;
        caxis([uncertainty_color_min uncertainty_color_max])
        title(h2, '\sigma_{\Delta x} (pix.)')
        title({'Moment of'; 'Correlation'})

        subplot(2,2,4)
        contourf(results.X, results.Y, sqrt(results.uncertainty2D.Ucsx.^2 + results.uncertainty2D.Ucsy.^2), uncertainty_contour_levels, 'edgecolor', 'none');
        colormap(flipud(gray))
        annotate_image(gcf, gca);
        h2 = colorbar;
        caxis([uncertainty_color_min uncertainty_color_max])
        title(h2, '\sigma_{\Delta x} (pix.)')
        title({'Correlation'; 'Statistics'})
        
        
        set(gcf, 'Position', [132    34   646   851]);

        if save_figures
            save_figure_to_png_eps_fig(figure_write_directory, 'prana-im-mc-contours', [true, false, false]);    
        end

    end
end

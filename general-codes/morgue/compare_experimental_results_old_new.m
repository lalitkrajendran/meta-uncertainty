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

% directory containing results for old processing
old_results_directory = fullfile('/scratch/shannon/a/aether/Projects/PIV_PTV_Uncertainty_Quantification/Planar_uncertainty_work/Processing/07_31_2018/');

% directory containing results for paper processing
paper_results_directory = fullfile('/scratch/shannon/a/aether/Projects/PIV_PTV_Uncertainty_Quantification/Planar_uncertainty_work/Results/07_30_2018/');

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
        
    for dataset_index = 5 %1:num_datasets
        % extract data set name
        dataset_name = dataset_name_array{dataset_index};
        fprintf('dataset name: %s\n', dataset_name);
        
        %% load new results
        
        % results directory for current data set
        current_results_directory = fullfile(top_write_directory, dataset_name, ['WS' num2str(window_size_index)]);

        % directory containing vectors
        vectors_directory = fullfile(current_results_directory, 'vectors');

        % get list of mat files containing results for the desired pass
        [files_new, num_files_new] = get_directory_listing(vectors_directory, ['*pass' num2str(pass_number) '*.mat']);
        
        % load sample result from new processing
        results_new = load(fullfile(files_new(1).folder, files_new(1).name));
        
        % calculate size of the new results
        [num_rows_new, num_cols_new] = size(results_new.X);
        
        %% load old results from processing
        
        % directory containing old vectors
        old_vectors_directory = fullfile(old_results_directory, ['WS' num2str(window_size_array(window_size_index))], dataset_name);

        % get list of mat files containing results for the desired pass
        [files_old, num_files_old] = get_directory_listing(old_vectors_directory, ['*pass' num2str(pass_number) '*.mat']);
        
        % load sample result from new processing
        results_old = load(fullfile(files_old(1).folder, files_old(1).name));

        % calculate size of the old results
        [num_rows_old, num_cols_old] = size(results_old.X);

        %% load old results from MC paper
        
        % directory containing old vectors
        paper_vectors_directory = fullfile(paper_results_directory, ['WS' num2str(window_size_array(window_size_index))], 'SOC_gradient', 'matfiles');
        
        % load sample result from new processing
        results_paper = load(fullfile(paper_vectors_directory, [dataset_name '.mat']));

        % calculate size of the old results
        [num_rows_paper, num_cols_paper] = size(results_paper.Xp);
        num_files_paper = size(results_paper.Up, 3);
                
        %% compare results
        
        % display size of new results
        fprintf('NEW: %d, %d, %d\n', num_rows_new, num_cols_new, num_files_new);
        
        % display size of old results
        fprintf('OLD: %d, %d, %d\n', num_rows_old, num_cols_old, num_files_old);

        % display size of paper results
        fprintf('PAPER: %d, %d, %d\n', num_rows_paper, num_cols_paper, num_files_paper);
        
    end
end
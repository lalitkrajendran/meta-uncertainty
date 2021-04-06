clear
close all
clc

restoredefaultpath;
addpath(genpath('/scratch/shannon/c/aether/Projects/BOS/general-codes/matlab-codes/'));
setup_default_settings;
% addpath ../prana-uncertainty-average-dc-new-im-2-save-deform-cs/
addpath ../prana/

%% processing settings

% directory containing images
top_image_directory = fullfile('/scratch/shannon/a/aether/Projects/PIV_PTV_Uncertainty_Quantification/Planar_uncertainty_work/Images/experiment/');

% array of data set names
dataset_name_array = {'PivChal03B'; 'PivChal05B'; 'stagnation_flow'; 'Vortex_Ring'; 'Jetdata'};
num_datasets = numel(dataset_name_array);

% directory containing sample job files
sample_job_file_directory = '/scratch/shannon/a/aether/Projects/PIV_PTV_Uncertainty_Quantification/uncertainty-combination/analysis/data/sample-job-files/';

% directory where results of this analysis are to be saved
% top_write_directory = fullfile('/scratch/shannon/a/aether/Projects/PIV_PTV_Uncertainty_Quantification/uncertainty-combination/analysis/results/planar-dataset/', 'experiment-new-grid-buffer');
top_write_directory = fullfile('/scratch/shannon/a/aether/Projects/PIV_PTV_Uncertainty_Quantification/uncertainty-combination/analysis/results/planar-dataset/', 'experiment-new');
if ~exist(top_write_directory, 'dir')
    mkdir(top_write_directory);
end

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
        mkdir_c(current_results_directory);
        
        % load job file
        sample_job_file = load(fullfile(current_results_directory, 'jobfile-lsg.mat'));
        
        % extract data structure containing job parameters
        Data = sample_job_file.Data;

        % directory to save vectors
        Data.outdirec = fullfile(current_results_directory, 'vectors-new-lsg');
        mkdir_c(Data.outdirec);

        % turn on ppr and mi uncertainty calc
        Data.PIV4.ppruncertainty = 1;
        Data.PIV4.miuncertainty = 1;

        % save jobfile
        save(fullfile(current_results_directory, 'jobfile-lsg.mat'), 'Data');
        
        % turn off parallel processing if just one image pair needs to be
        % processed
        if str2double(Data.imfend) == 1
            Data.par = '0';
        end
        
        % call prana
        pranaPIVcode(Data);        
    end
end

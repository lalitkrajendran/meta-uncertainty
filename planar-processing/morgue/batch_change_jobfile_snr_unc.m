clear
close all
clc

restoredefaultpath;
addpath(genpath('/scratch/shannon/c/aether/Projects/BOS/general-codes/matlab-codes/'));
setup_default_settings;
addpath ../prana-uncertainty-average-dc-new-im-2-save-deform-cs/
% addpath ../prana/
%% processing settings

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

%% run processing

for window_size_index = 1:2
    fprintf('WS: %d\n', window_size_index);
        
    for dataset_index = 1:num_datasets
        % extract data set name
        dataset_name = dataset_name_array{dataset_index};
        fprintf('dataset name: %s\n', dataset_name);
        
        % results directory for current data set
        current_results_directory = fullfile(top_write_directory, dataset_name, ['WS' num2str(window_size_index)]);
        
        % load job file
        sample_job_file = load(fullfile(current_results_directory, 'jobfile.mat'));
        
        % extract data structure containing job parameters
        Data = sample_job_file.Data;
        
        %%
        field_names = fieldnames(Data);
        for field_name_index = 1:numel(field_names)
            % extract field name
            field_name = field_names{field_name_index};
            % check if the field name corresponds to a pass and copy over
            if contains(field_name, 'PIV') && ~contains(field_name, 'run') && ~contains(field_name, 'error')
                fprintf('%s\n', field_name);
                % calculate number of passes
                num_passes = str2double(Data.passes);
                if contains(field_name, num2str(num_passes-1))
                    % turn on uncertainty estimate on the third pass
                    Data.(field_name).uncertaintyestimate = '1';
                    Data.(field_name).imuncertainty = '1';
                    Data.(field_name).mcuncertainty = '1';
                    Data.(field_name).csuncertainty = '1';
                end
            end
        end
        
        %% make other changes
        % add option to save deform images
        Data.SaveIMdeform = '1';
        
        % add option save cs uncertainty on last pass
        Data.PIV4.csuncertainty = '1';
        
        %% save
        % make new results directory
        new_results_directory = fullfile(top_write_directory, dataset_name, ['WS' num2str(window_size_index)], 'new-grid-buffer');
        mkdir_c(new_results_directory);
        % save jobfile
        save(fullfile(new_results_directory, 'jobfile.mat'), 'Data');
                
    end
end

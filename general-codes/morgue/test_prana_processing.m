clear
close all
clc

restoredefaultpath;
% addpath prana-uncertainty-average-dc-new-im-2-save-deform-cs/
addpath ../prana/

% load sample job file
job_file_name = '/scratch/shannon/a/aether/Projects/PIV_PTV_Uncertainty_Quantification/uncertainty-combination/analysis/results/prana-test-cs-buffer=12/jobfile.mat';
sample_job_file = load(job_file_name);

% extract processing setting structure
Data = sample_job_file.Data;

% adjust output directory
Data.outdirec = '/scratch/shannon/a/aether/Projects/PIV_PTV_Uncertainty_Quantification/uncertainty-combination/analysis/results/prana-test-cs-buffer=12/';
% adjust grid buffer
Data.PIV1.gridbuf = '12,12';
Data.PIV2.gridbuf = '12,12';
Data.PIV3.gridbuf = '12,12';
Data.PIV4.gridbuf = '12,12';

% turn on cs uncertainty for last pass
Data.PIV1.imuncertainty = '0';
Data.PIV2.imuncertainty = '0';
Data.PIV3.imuncertainty = '0';
Data.PIV4.imuncertainty = '0';

Data.PIV1.mcuncertainty = '0';
Data.PIV2.mcuncertainty = '0';
Data.PIV3.mcuncertainty = '0';
Data.PIV4.mcuncertainty = '0';

Data.PIV1.csuncertainty = '0';
Data.PIV2.csuncertainty = '0';
Data.PIV3.csuncertainty = '0';
Data.PIV4.csuncertainty = '1';

% turn on saving deform images
Data.SaveIMdeform = '1';
% 
% % save jobfile
% save(job_file_name, 'Data');

% run prana
pranaPIVcode(Data);




% script to compare mean flows from prana and davis processing 
% for the stereo vortex ring dataset

clear
close all
clc

restoredefaultpath;
addpath(genpath('/scratch/shannon/c/aether/Projects/BOS/general-codes/matlab-codes/'));
addpath('../prana-uncertainty-average-dc-new-im-2-save-deform-cs/');
addpath('../stereo_uncertainty_codes_packaged/');
setup_default_settings;

% ============================
% directory settings
% ============================
% top level read directory
top_read_directory = '/scratch/shannon/a/aether/Projects/PIV_PTV_Uncertainty_Quantification/Stereo_uncertainty_work/';
% top level write directory
top_write_directory = '/scratch/shannon/a/aether/Projects/PIV_PTV_Uncertainty_Quantification/uncertainty-combination/analysis/results/stereo-dataset/cam13/With_selfcal/';
% directory containing old results
code_package_directory = '../stereo_uncertainty_codes_packaged/';

% ============================
% processing settings
% ============================
% stereo reconstruction type
rectype = 'Willert';
% pass number of results
pass_index = 4;

% ============================
% plot settings
% ============================
% file index to compare
file_index = 1;
% column index for the line plot
prana_col = 55;
davis_col = 55;

% displacement scale factor
quiver_scale_factor = 10.0;
% number of vectors to skip
vector_skip = 2;
% quiver line width
quiver_width = 1.5;
% save figure? (true/false)
save_figures = 0;
% directory to save figures
figure_save_directory = fullfile(top_write_directory, 'figures');
mkdir_c(figure_save_directory);

% ============================
% Load Stereo calibration job
% ============================
stereojob = load(fullfile(code_package_directory,'Prana_stereo_job_after_selfcal.mat'));
stdjob = stereojob.stdjob;

% extract Calibration job
caljob = stdjob.caljobfile.datasave.caljob;
calmat = [caljob.aXcam1 caljob.aYcam1 caljob.aXcam2 caljob.aYcam2];

% ============================
% load prana solution
% ============================
% directory containing processed result
vectors_directory = fullfile(top_write_directory, rectype, ['Camera',num2str(caljob.camnumber(1)),'Camera',num2str(caljob.camnumber(2)),'_3Cfields',filesep]);
% get list of vector fields for the specified pass
[vector_files, num_vector_files] = get_directory_listing(vectors_directory, ['piv*pass_' num2str(pass_index, '%d') '*.mat']);
% load mat files
[vector_results, num_vector_results] = load_directory_data(vectors_directory, ['piv*pass_' num2str(pass_index, '%d') '*.mat']);

% ============================
% load davis solution
% ============================
Davissol = load(fullfile(code_package_directory,'Davis_processed_result.mat'));
Davis13 = Davissol.Davis13;

% ============================
% plot contours side by side
% ============================
figure
% quiver
subplot(1, 3, 1)
quiver_skip(gcf, gca, vector_results(file_index).X*1e3, vector_results(file_index).Y*1e3, vector_results(file_index).U, ...
        vector_results(file_index).V, vector_skip, quiver_scale_factor, quiver_width);
hold on
quiver_skip(gcf, gca, Davis13.Xw*1e3, Davis13.Yw*1e3, Davis13.Uw(:, :, file_index), ...
        Davis13.Vw(:, :, file_index), vector_skip, quiver_scale_factor, quiver_width);
set_axes(gca);
legend('Prana', 'Davis', 'location', 'northoutside')
% xlabel('X (mm)');
% ylabel('Y (mm)');

% line plot of u
subplot(1, 3, 2)
plot(vector_results(file_index).U(:, prana_col), vector_results(file_index).Y(:, 1)*1e3, 'o')
hold on
plot(Davis13.Uw(:, davis_col, file_index), Davis13.Yw(:, 1)*1e3, 'o')
xlabel('U')
ylabel('Y')
xlim([-0.1 0.2])

% line plot of V
subplot(1, 3, 3)
plot(vector_results(file_index).V(:, prana_col), vector_results(file_index).Y(:, 1)*1e3, 'o')
hold on
plot(Davis13.Vw(:, davis_col, file_index), Davis13.Yw(:, 1)*1e3, 'o')
xlabel('V')
ylabel('Y')
xlim([-0.1 0.2])

set(gcf, 'resize', 'off');
set(gcf, 'Position', [65          41        1167         660]);
% quiver(vector_results(1).X*1e3, vector_results(1).Y*1e3, vector_results(1).U * quiver_scale_factor, vector_results(1).V * quiver_scale_factor, 'autoscale', 'off');

% ============================
% make line plot of a slice through the vortex ring
% ============================

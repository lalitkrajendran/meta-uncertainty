% This script propagates the planar uncertainties obtained from 
% jack-knife resampling and calculates the 3c uncertainties.
% It then uses the distribution of the 3c uncertainties to calculate
% weights and combined stereo uncertainties

%%
clear
close all
clc

restoredefaultpath;
addpath(genpath('/scratch/shannon/c/aether/Projects/BOS/general-codes/matlab-codes/'));
addpath('../prana-uncertainty-average-dc-new-im-2-save-deform-cs/');
addpath('../general-codes/');
addpath('../stereo_uncertainty_codes_packaged/');
addpath('../CompPD/')

setup_default_settings;
% dbstop if error

% ============================
%% experiment settings
% ============================
% seconds per frame
spf = 0.001;
% No. of frames
num_frames = 50;   
% Starting frame
fstart = 24;
% camera numbers
camera_numbers = [2, 4];
num_cameras = numel(camera_numbers);

% ====================================
%% read/write settings
% ====================================
% top level read directory
top_read_directory = '/scratch/shannon/a/aether/Projects/PIV_PTV_Uncertainty_Quantification/Stereo_uncertainty_work/';
% top level write directory
% top_write_directory = '/scratch/shannon/a/aether/Projects/PIV_PTV_Uncertainty_Quantification/uncertainty-combination/analysis/results/stereo-dataset/cam13/With_selfcal_old/';
% top_write_directory = '/scratch/shannon/a/aether/Projects/PIV_PTV_Uncertainty_Quantification/uncertainty-combination/analysis/results/stereo-dataset/cam13/With_selfcal/';
top_write_directory = fullfile('/scratch/shannon/a/aether/Projects/PIV_PTV_Uncertainty_Quantification/uncertainty-combination/', '/analysis/results/stereo-dataset/', ['cam' num2str(camera_numbers(1)) num2str(camera_numbers(2))], '/With_selfcal/');
mkdir_c(top_write_directory);
% Load 2d job file 
% job_settings = load(fullfile(top_write_directory, 'cam13_VR_prana_fulljob_withselfcal.mat'));
job_settings = load(fullfile(top_write_directory, ['cam' num2str(camera_numbers(1)) num2str(camera_numbers(2)) '_VR_prana_fulljob_withselfcal.mat']));
% extract calibration job file
caljobfile = job_settings.caljobfile;
% extract 2d job file
planarjob = job_settings.planarjob;

% ====================================
%% processing settings
% ====================================
% % camera numbers
% camera_numbers = [1, 3];
% num_cameras = numel(camera_numbers);
% stereo reconstruction type
rectype = 'Willert';
% pass number
pass_number = 4;
% array of uncertainty methods
individual_method_array = {'IM'; 'MC'; 'CS'};
num_individual_methods = numel(individual_method_array);
% name of the uncertainty methods in prana
uncertainty_method_prana_names = {'Uim'; 'MC'; 'Ucs'};
% start frame
start_frame = 25;
% component names
component_names = {'x'; 'y'; 'z'};
num_components = numel(component_names);

% ============================
%% resampling settings
% ============================
% resampling method names
resampling_method_name_array = {'remove-paired'; 'remove-random'; 'add-random'};
% resampling_method_name_array = {'add-random'};
num_resampling_methods = numel(resampling_method_name_array);
% number of particles to remove
percentage_particles_remove_array = 0.05:0.05:0.25; %0.025:0.025:0.25; 
num_ppr = numel(percentage_particles_remove_array);
% number of resampling trials
num_resampling_trials = 1e2;
% num overall trials
num_trials = 1e3;

% ========================
%% weight settings
% ========================
% number of bins for weights
num_bins_w = 30;
% bins for weights
bins_w = linspace(0, 1, num_bins_w);
% metric type ('rms', 'iqr', 'std')
metric_name = 'iqr';

% ====================================
% combination method settings
% ====================================
% array of uncertainty combination model to be used
% 'unweighted', 'global-weight-var', 'global-weight-std',
% 'local-weight-var', 'local-weight-std'
% combination_methods = {'unwt'; 'var'; 'pd'; 'entropy'}; %'prob'};
combination_methods = {'unwt'; 'pd-var'; 'entropy'}; %'prob'};
num_combination_methods = numel(combination_methods);

% ===============================
%% analysis settings
% ===============================
% minimum allowable error (pix.)
error_min = 1e-4;
% maximum allowable error (pix.)
error_max = [0.5, 0.5, 1.5];
% number of bins for histogram
num_bins = 30;
% bins for histograms
bins = cell(1, 4);
for i = 1:3
    bins{i} = linspace(error_min, error_max(i), num_bins);
end

% ====================================
%% load planar results
% ====================================
fprintf('Loading planar results\n');

results_planar_all = cell(1, 2);
jobfile_all = cell(1, 2);
files_im1 = cell(1, 2);
files_im2 = cell(1, 2);

% loop through cameras
for camera_index = 1:num_cameras
    fprintf('Camera: %d\n', camera_index);
    % ====================================
    %% Load data
    % ====================================
    % results directory for the current camera
    current_results_directory = fullfile(top_write_directory, rectype, ['Camera', num2str(camera_index), filesep], 'vectors');
    % load results for vectors, errors and uncertainties
    [results_planar_all{camera_index}, num_snapshots] = load_directory_data(current_results_directory, ['VR*pass' num2str(pass_number, '%d') '*.mat']);
    % extract number of grid point
    [num_rows, num_cols] = size(results_planar_all{camera_index}{1}.X);
    num_grid_points = num_rows * num_cols;    
end

% ====================================
% extract job file properties
% ====================================
% extract job file
jobfile = job_settings.job1;

% ====================================
%% load stereo results
% ====================================
fprintf('Loading stereo results\n');
% directory containing stereo results
stereo_results_directory = fullfile(top_write_directory, rectype, 'Camera1Camera2_3Cfields',filesep);
% load results for vectors, errors and uncertainties
[results_stereo_all, num_snapshots] = load_directory_data(fullfile(stereo_results_directory, 'vectors'), ['piv*pass_' num2str(pass_number, '%d') '*.mat']);
    
% load errors
fprintf('Loading all errors into memory\n');
filename = fullfile(stereo_results_directory, 'errors.mat');
errors_all = load(filename);

% load uncertainties
fprintf('Loading all uncertainties into memory\n');
filename = fullfile(stereo_results_directory, 'uncertainties.mat');
uncertainties_all = load(filename);

% ============================
% Extract calibration details
% ============================
% Magnification in mm/pix
scaling = job_settings.scaling.wil;
mx = scaling.xscale;
my = scaling.yscale;
% Here the results are for camera1camera3 pair
cam1 = '1';
cam2 = '3';
if strcmp(cam1,'1') && strcmp(cam2,'3')
    mz = scaling.yscale;
elseif strcmp(cam1,'2') && strcmp(cam2,'4')
    mz = scaling.xscale;
end

% ===============================
%% load resampling results
% ===============================

% directory to load results for this case
current_read_directory = fullfile(stereo_results_directory, 'resampling-new', ...
                                ['trials=' num2str(num_trials, '%d'), '_nrs' num2str(num_resampling_trials) '-prewindow']);

% load monte carlo results
load(fullfile(current_read_directory, 'monte_carlo_results.mat'));

% ===============================
%% propagate resampling results and calculate stereo uncertainty
% ===============================
unc_stereo_individual = cell(1, num_trials);
unc_stereo_resampled = cell(num_trials, num_ppr, num_resampling_methods);
unc_stereo_combined = cell(1, num_trials);
pdf_unc_resampled = cell(1, num_trials);
weights = cell(1, num_trials);
unc_ratio = cell(1, num_trials);
unc_ratio_fit = cell(1, num_trials);

for trial_index = 1:num_trials
    fprintf('trial: %d\n', trial_index);
    % ================================================
    % extract velocities, errors and uncertainties
    % ================================================
    % extract snapshot for this trial
    f_trial = f_trial_all(trial_index);

    % extract grid points for this trial
    r_trial = r_trial_all(trial_index);
    c_trial = c_trial_all(trial_index);

    % extract results
    results_planar = results_planar_all{1}{f_trial};    

    % extract 2c velocities    
    U1 = results_planar_all{1}{f_trial}.U(r_trial, c_trial);
    V1 = results_planar_all{1}{f_trial}.V(r_trial, c_trial);
    U2 = results_planar_all{2}{f_trial}.U(r_trial, c_trial);
    V2 = results_planar_all{2}{f_trial}.V(r_trial, c_trial);
    
    % extract 3c velocities
    w_current = (spf/mx) * 1e3 * results_stereo_all{f_trial}.W(r_trial, c_trial);

    % extract calibration angles
    tanalpha1 = uncertainties_all.tanalpha1(r_trial, c_trial);
    tanalpha2 = uncertainties_all.tanalpha2(r_trial, c_trial);
    tanbeta1 = uncertainties_all.tanbeta1(r_trial, c_trial);
    tanbeta2 = uncertainties_all.tanbeta2(r_trial, c_trial);

    % extract angle uncertainties
    Un_alpha1 = uncertainties_all.Un_alpha1(r_trial, c_trial);
    Un_alpha2 = uncertainties_all.Un_alpha2(r_trial, c_trial);
    Un_beta1 = uncertainties_all.Un_beta1(r_trial, c_trial);
    Un_beta2 = uncertainties_all.Un_beta2(r_trial, c_trial);
    
    % ================================================
    % calculate stereo uncertainties from ORIGINAL planar uncertainties
    % ================================================
    unc_stereo_individual{trial_index} = cell(1, num_individual_methods);
    for individual_method_index = 1:num_individual_methods
        % ================================================
        % extract planar uncertainties from resampling
        % ================================================
        % camera 1
        Unu1 = unc_planar_individual{1, trial_index}{individual_method_index}.x;
        Unv1 = unc_planar_individual{1, trial_index}{individual_method_index}.y;
        % camera 2
        Unu2 = unc_planar_individual{2, trial_index}{individual_method_index}.x;
        Unv2 = unc_planar_individual{2, trial_index}{individual_method_index}.y;

        % ================================================
        % propagate uncertainties
        % ================================================
        [Un_u1, Un_v1, Un_w1] = stereo_uncertainty_propagation(Unu1, Unv1, Unu2, Unv2, ...
                                                                U1, V1, U1, V1, w_current, ...
                                                                Un_alpha1, Un_alpha2, Un_beta1, Un_beta2, ...
                                                                tanalpha1, tanalpha2, tanbeta1, tanbeta2, ...
                                                                mx, my);
        % store uncertainties
        unc_stereo_individual{trial_index}{individual_method_index}.x = Un_u1;
        unc_stereo_individual{trial_index}{individual_method_index}.y = Un_v1;
        unc_stereo_individual{trial_index}{individual_method_index}.z = Un_w1;
    end    
    
    % ================================================
    % calculate stereo uncertainties from RESAMPLED planar uncertainties
    % ================================================
    % ==========================
    % loop through resampling methods
    % ==========================
    weights{trial_index} = cell(1, num_resampling_methods);
    unc_ratio{trial_index} = cell(1, num_resampling_methods);
    unc_ratio_fit{trial_index} = cell(1, num_resampling_methods);
    unc_stereo_combined{trial_index} = cell(1, num_resampling_methods);

    for resampling_method_index = 1:num_resampling_methods
        % fprintf('resampling method index: %d\n', resampling_method_index);
        resampling_method_name = resampling_method_name_array{resampling_method_index};        
        % % copy results
        % unc_planar_resampling{camera_index, trial_index, particle_remove_index, resampling_method_index} = unc;
        % snr_planar_resampling{camera_index, trial_index, particle_remove_index, resampling_method_index} = snr;
        % ================================================
        %% loop through particle removal percentages
        % ================================================
        for particle_remove_index = 1:num_ppr
            % current percentage to remove
            percentage_particles_remove = percentage_particles_remove_array(particle_remove_index);
            unc_stereo_resampled{trial_index, particle_remove_index, resampling_method_index} = cell(1, num_individual_methods);
            % fprintf('particle remove percentage: %.2f\n', percentage_particles_remove);                    
            % ================================================
            %% loop through individual methods
            % ================================================
            for individual_method_index = 1:num_individual_methods
                % initialize
                unc_stereo_resampled{trial_index, particle_remove_index, resampling_method_index}{individual_method_index} = struct;
                unc_stereo_resampled{trial_index, particle_remove_index, resampling_method_index}{individual_method_index}.x = nans(1, num_resampling_trials);
                unc_stereo_resampled{trial_index, particle_remove_index, resampling_method_index}{individual_method_index}.y = nans(1, num_resampling_trials);
                unc_stereo_resampled{trial_index, particle_remove_index, resampling_method_index}{individual_method_index}.z = nans(1, num_resampling_trials);
                % ================================================
                % loop through resampling trials 
                % ================================================
                for resampling_trial_index = 1:num_resampling_trials                    
                    % ================================================
                    % extract planar uncertainties from resampling
                    % ================================================
                    % camera 1
                    Unu1 = unc_planar_resampling{1, trial_index, particle_remove_index, resampling_method_index}.([lower(individual_method_array{individual_method_index}) 'x'])(resampling_trial_index);
                    Unv1 = unc_planar_resampling{1, trial_index, particle_remove_index, resampling_method_index}.([lower(individual_method_array{individual_method_index}) 'y'])(resampling_trial_index);
                    % camera 2
                    Unu2 = unc_planar_resampling{2, trial_index, particle_remove_index, resampling_method_index}.([lower(individual_method_array{individual_method_index}) 'x'])(resampling_trial_index);
                    Unv2 = unc_planar_resampling{2, trial_index, particle_remove_index, resampling_method_index}.([lower(individual_method_array{individual_method_index}) 'y'])(resampling_trial_index);

                    % ================================================
                    % propagate uncertainties
                    % ================================================
                    [Un_u1, Un_v1, Un_w1] = stereo_uncertainty_propagation(Unu1, Unv1, Unu2, Unv2, ...
                                                                            U1, V1, U2, V2, w_current, ...
                                                                            Un_alpha1, Un_alpha2, Un_beta1, Un_beta2, ...
                                                                            tanalpha1, tanalpha2, tanbeta1, tanbeta2, ...
                                                                            mx, my);
                    % store uncertainties
                    unc_stereo_resampled{trial_index, particle_remove_index, resampling_method_index}{individual_method_index}.x(resampling_trial_index) = Un_u1;
                    unc_stereo_resampled{trial_index, particle_remove_index, resampling_method_index}{individual_method_index}.y(resampling_trial_index) = Un_v1;
                    unc_stereo_resampled{trial_index, particle_remove_index, resampling_method_index}{individual_method_index}.z(resampling_trial_index) = Un_w1;
                end
            end
        end    

        % ============================
        % calculate weights
        % ============================
        % [wt_trials{resampling_method_index, trial_index}, unc_ratio{resampling_method_index, trial_index}, unc_ratio_fit{resampling_method_index, trial_index}] = calculate_weights_from_resampled_uncertainty_ratios(unc_sub_trials{trial_index}, {unc_resampling_trials{trial_index, :, resampling_method_index}}, ...
        % metric_name, individual_method_array, percentage_particles_remove_array, components);

        [weights{trial_index}{resampling_method_index}, unc_ratio{trial_index}{resampling_method_index}, unc_ratio_fit{trial_index}{resampling_method_index}] = calculate_weights_from_resampled_uncertainty_ratios(unc_stereo_individual{trial_index}, {unc_stereo_resampled{trial_index, :, resampling_method_index}}, ...
                                                                    metric_name, individual_method_array, percentage_particles_remove_array, component_names, num_resampling_trials, 'cell');

        % ============================
        % calculate combined stereo uncertainty
        % ============================
        unc_stereo_combined{trial_index}{resampling_method_index} = struct;
        for component_index = 1:num_components
            component_name = component_names{component_index};
            unc_indiv = [unc_stereo_individual{trial_index}{1}.(component_name), unc_stereo_individual{trial_index}{2}.(component_name), ...
                            unc_stereo_individual{trial_index}{3}.(component_name)];

            unc_stereo_combined{trial_index}{resampling_method_index}.(component_name) = unc_indiv * weights{trial_index}{resampling_method_index}.(component_name)';
        end
    end     
end

shut_down_parpool();
% ==================
% save results
% ==================
fprintf('saving results to file\n');
filename = fullfile(current_read_directory, 'combined_uncertainties_3c.mat');
save(filename, 'weights', 'unc_ratio', 'unc_ratio_fit', ... 
                'unc_stereo_individual', 'unc_stereo_resampled', 'unc_stereo_combined');


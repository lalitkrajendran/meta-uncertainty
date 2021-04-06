% This script performs resampling calculations for a range of particle removal
% percentage and plots the change in uncertainties for different methods

%%
clear
close all
clc

restoredefaultpath;
addpath(genpath('/scratch/shannon/c/aether/Projects/BOS/general-codes/matlab-codes/'));
addpath(genpath('/scratch/shannon/c/aether/Projects/BOS/general-codes/piv-image-generation/'));
% addpath('../prana-uncertainty-average-dc-new-im-2-save-deform-cs/');
addpath ../prana/
addpath ../general-codes/

setup_default_settings;

% ============================
%% read/write settings
% ============================
% window resolution
window_resolution_array = [64, 32];
num_window_resolution = numel(window_resolution_array);
% pass number
pass_number = 4;
% directory containing images
top_image_directory = fullfile('/scratch/shannon/a/aether/Projects/PIV_PTV_Uncertainty_Quantification/Planar_uncertainty_work/Images/');
% directory containing files to be read
top_read_directory = fullfile('/scratch/shannon/a/aether/Projects/PIV_PTV_Uncertainty_Quantification/uncertainty-combination/analysis/results/planar-dataset', 'experiment-new');
% array of data set names
dataset_name_array = {'PivChal03B'; 'PivChal05B'; 'stagnation_flow'; 'Vortex_Ring'; 'Jetdata'};
num_datasets = numel(dataset_name_array);
% array of uncertainty methods
individual_method_array = {'IM'; 'MC'; 'CS'};
num_individual_methods = numel(individual_method_array);
% array of snr metrics
snr_metric_array = {'PPR'; 'MI'};
num_snr_methods = numel(snr_metric_array);
% directory where results of this analysis are to be saved
top_write_directory = fullfile('/scratch/shannon/a/aether/Projects/PIV_PTV_Uncertainty_Quantification/uncertainty-combination/analysis/results/planar-dataset', 'monte-carlo', 'individual-datasets');
mkdir_c(top_write_directory);
% number of trials
num_trials = 1; %000; %1; %1e1;

% ============================
%% resampling settings
% ============================
% number of particles to remove
% percentage_particles_remove_array = 0:0.05:0.25;
percentage_particles_remove_array = 0.05:0.05:0.25;
num_ppr = numel(percentage_particles_remove_array);
% number of resampling trials
num_resampling_trials = 1e3;
% display_resampled_images? (true/false)
display_resampled_images = 0;
% % case name
% resampling_case_name = ['ppr' num2str(percentage_particles_remove, '%.2f') '_nrs' num2str(num_resampling_trials) '-prewindow-randpart-snr-newproc'];
% resampling method names
resampling_method_names = {'Removing Paired Particles'; 'Removing Random Particles'; 'Adding Random Particles'};
num_resampling_methods = numel(resampling_method_names);

% ============================
%% particle identification settings
% ============================
% intensity threshold for particle identification
intensity_threshold = 10;
% particle diameter (pix.)
d_p = 3;
% particle sizing settings
sizeprops = struct;
sizeprops.method = 'IWC';
sizeprops.p_area = 2; % this will be modified for each image
sizeprops.sigma = 4;
sizeprops.errors = 0;
sizeprops.save_dir = [];
% display id results?
display_id_results = 0;

% ========================
%% plot settings
% ========================
% display figures?
display_figures = 1; %0;
% user screen resolution
user_screen_resolution = 113;
% save figure? (true/false)
save_figures = 0; %1;
% symbols
symbols = {'o', '^', 'v'};
% colors
colors = lines(3);
% line symbols
line_symbols = {'-'; ':'; '-.'};
% range of displacements to be displayed in the contour plots
displacement_color_min = [0, 0, -0.25, -2, 0];
displacement_color_max = [15, 5, 0.25, 2, 5]; 
% displacement_contour_levels = linspace(displacement_color_min, displacement_color_max, 100);

% ============================
%% pre-load all dataset
% ============================
fprintf('Loading all datasets into memory\n');

results_all = cell(num_window_resolution, num_datasets);
jobfile_all = cell(num_window_resolution, num_datasets);
files_im1 = cell(num_window_resolution, num_datasets);
files_im2 = cell(num_window_resolution, num_datasets);

% loop through window resolutions
for window_resolution_index = 1:num_window_resolution
    % loop through datasets
    for dataset_index = 1:num_datasets
        % name of the data set
        dataset_name = dataset_name_array{dataset_index};        
        fprintf('Dataset: %s\n', dataset_name);

        % ============================
        %% load data
        % ============================
        % results directory for current data set
        current_results_directory = fullfile(top_read_directory, dataset_name, ['WS' num2str(window_resolution_index)]);

        % directory containing vectors
        vectors_directory = fullfile(current_results_directory, 'vectors-new');

        % load results for vectors, errors and uncertainties
        results_all{window_resolution_index, dataset_index} = load_directory_data(vectors_directory, ['*pass' num2str(pass_number) '*.mat']);

        % load jobfile
        jobfile_all{window_resolution_index, dataset_index} = load(fullfile(current_results_directory, 'jobfile.mat'));
        
        % ============================
        %% load listing of deformed images
        % ============================        
        % directory containing deformed images for current data set
        deformed_images_directory = fullfile(vectors_directory, 'imDeform');
        
        % get list of im1 files in the directory
        files_im1{window_resolution_index, dataset_index} = get_directory_listing(deformed_images_directory, 'PIV*im1d*.mat');
        % get list of im2 files in the directory
        files_im2{window_resolution_index, dataset_index} = get_directory_listing(deformed_images_directory, 'PIV*im2d*.mat');        
    end
end

% ============================
%% load errors for all datasets
% ============================
fprintf('Loading all errors into memory\n');
errors_all = cell(num_window_resolution, num_datasets);

% loop through window resolutions
for window_resolution_index = 1:num_window_resolution
    % loop through datasets
    for dataset_index = 1:num_datasets
        % name of the data set        
        dataset_name = dataset_name_array{dataset_index};
        
        fprintf('Dataset: %s\n', dataset_name);
        
        %% Load data

        % results directory for current data set
        current_results_directory = fullfile(top_read_directory, dataset_name, ['WS' num2str(window_resolution_index)]);

        % load results for vectors, errors and uncertainties
        errors_all{window_resolution_index, dataset_index} = load(fullfile(current_results_directory, 'errors-new.mat'));        
    end
end

fprintf('running monte-carlo\n');
% set seed for random number generator
rng(0);

X_all = cell(num_window_resolution, num_datasets);
Y_all = cell(num_window_resolution, num_datasets);
unc_mean_all = cell(num_window_resolution, num_datasets);
d_unc_mean_all = cell(num_window_resolution, num_datasets);
snr_mean_all = cell(num_window_resolution, num_datasets);

% start timer
tic
% ============================
%% loop through window resolutions
% ============================
for window_resolution_index = 1:num_window_resolution
    fprintf('WS: %d\n', window_resolution_index);

    % ============================
    %% loop through datasets
    % ============================
    for dataset_index = 1:num_datasets
        % dataset name
        dataset_name = dataset_name_array{dataset_index};
        fprintf('dataset: %s\n', dataset_name);

        % ================================================
        % load results
        % ================================================
        % directory to save results for this case
        % current_write_directory = fullfile(top_write_directory, dataset_name, ['WS' num2str(window_resolution_index)], 'resampling-study');
        current_write_directory = '/scratch/shannon/a/aether/Projects/PIV_PTV_Uncertainty_Quantification/uncertainty-combination/analysis/results/planar-dataset/monte-carlo/individual-datasets/PivChal03B/WS1/weights-rms-change-study';
        mkdir_c(current_write_directory);
        % file name
        filename = ['resampling-statistics-nt=' num2str(num_trials, '%d') 'nr='  num2str(num_resampling_trials, '%d') '.mat'];
        % save results to file
        load(fullfile(current_write_directory, filename));
        
        % create directory to save figures
        if save_figures
            figure_save_directory = fullfile(current_write_directory, 'figures');
            mkdir_c(figure_save_directory);
        end

        % ================================================
        %% loop through trials
        % ================================================
        for trial_index = 1 %1:num_trials
            fprintf('trial_index: %d\n', trial_index);

            % ==========================
            % loop through resampling methods
            % ==========================
            for resampling_method_index = 1:num_resampling_methods
                fprintf('resampling method index: %d\n', resampling_method_index);            
                
                % ==========================
                % loop through uncertainty methods
                % ==========================
                for individual_method_index = 1:num_individual_methods
                    % method name
                    method_name = lower(individual_method_array{individual_method_index});
                    % extract original uncertainties
                    unc_sub = unc_sub_trials{trial_index}{individual_method_index}.x;
                    unc_ratio = nans(num_resampling_trials, num_ppr);
                    ppr_name_array = cell(1, num_ppr);
                    % ================================================
                    %% loop through particle removal percentages
                    % ================================================
                    for particle_remove_index = 1:num_ppr
                        % current percentage to remove
                        percentage_particles_remove = percentage_particles_remove_array(particle_remove_index);
                        fprintf('particle remove percentage: %.2f\n', percentage_particles_remove);
        
                        % extract resampled uncertainties
                        unc_r = unc_resampling_trials{trial_index, particle_remove_index, resampling_method_index}.([method_name 'x']);

                        % calculate ratio
                        unc_ratio(:, particle_remove_index) = unc_r ./ unc_sub;

                        % name for this case
                        ppr_name_array{particle_remove_index} = num2str(round(percentage_particles_remove * 100), '%d');
                    end

                    % make violin plot of ratio
                    subplot_index = (individual_method_index - 1) * num_resampling_methods + resampling_method_index;
                    figure(1)
                    subplot(num_individual_methods, num_resampling_methods, subplot_index)
                    % violinplot_single_symm(unc_ratio, particle_remove_index, 0.3, lines(1));
                    violins = violinplot(unc_ratio, ppr_name_array, 'violincolor', colors(1, :), 'showdata', false, 'shownotches', true, 'showmean', true, ...
                                        'showrms', true, 'edgecolor', [1, 1, 1], 'ViolinAlpha', 0.5);

                    % --------------------------
                    % annotate figure
                    % --------------------------                    
                    xlim([1 num_ppr+0.5])
                    ylim([0 2])

                    if individual_method_index == num_individual_methods
                        if resampling_method_index <= 2
                            xlabel('Particle Removal %');
                        else
                            xlabel('Particle Addition %');
                        end
                    else
                        set(gca, 'xticklabel', []);
                    end

                    if resampling_method_index == 1
                        ylabel('Uncertainty Ratio')                    
                    else
                        set(gca, 'yticklabel', []);
                    end

                    if individual_method_index == 1
                        title(resampling_method_names{resampling_method_index})
                    end                    

                    if resampling_method_index == 1
                        annotation('textbox', [0.02, 0.85 - 0.3 * (individual_method_index- 1), 0, 0], 'string', upper(method_name), 'fontsize', 18, 'fontweight', 'bold')
                    end
                end

                % ==========================
                % loop through snr methods
                % ==========================
                for snr_method_index = 1:num_snr_methods
                    % method name
                    method_name = snr_metric_array{snr_method_index};
                    % extract original uncertainties
                    snr_sub = snr_sub_trials{trial_index}.(method_name);
                    snr_ratio = nans(num_resampling_trials, num_ppr);
                    ppr_name_array = cell(1, num_ppr);

                    % ================================================
                    %% loop through particle removal percentages
                    % ================================================
                    for particle_remove_index = 1:num_ppr
                        % current percentage to remove
                        percentage_particles_remove = percentage_particles_remove_array(particle_remove_index);
                        fprintf('particle remove percentage: %.2f\n', percentage_particles_remove);
        
                        % extract resampled snr metric
                        snr_r = snr_resampling_trials{trial_index, particle_remove_index, resampling_method_index}.(method_name);

                        % calculate ratio
                        snr_ratio(:, particle_remove_index) = snr_r ./ snr_sub;

                        % name for this case
                        ppr_name_array{particle_remove_index} = num2str(round(percentage_particles_remove * 100), '%d');
                    end

                    % make violin plot of ratio
                    subplot_index = (snr_method_index - 1) * num_resampling_methods + resampling_method_index;
                    figure(2)
                    subplot(num_snr_methods, num_resampling_methods, subplot_index)
                    % violinplot_single_symm(unc_ratio, particle_remove_index, 0.3, lines(1));
                    violins = violinplot(snr_ratio, ppr_name_array, 'violincolor', colors(1, :), 'showdata', false, 'shownotches', true, 'showmean', true, ...
                                        'showrms', true, 'edgecolor', [1, 1, 1], 'ViolinAlpha', 0.5);
                                            
                    % --------------------------
                    % annotate figure
                    % --------------------------                    
                    xlim([1 num_ppr+0.5])
                    ylim([0 2])

                    if snr_method_index == num_snr_methods
                        if resampling_method_index <= 2
                            xlabel('Particle Removal %');
                        else
                            xlabel('Particle Addition %');
                        end
                    else
                        set(gca, 'xticklabel', []);
                    end

                    if resampling_method_index == 1
                        ylabel('SNR Ratio')                    
                    else
                        set(gca, 'yticklabel', []);
                    end

                    if snr_method_index == 1
                        title(resampling_method_names{resampling_method_index})
                    end                    

                    if resampling_method_index == 1
                        annotation('textbox', [0.02, 0.8 - 0.5 * (snr_method_index - 1), 0, 0], 'string', upper(method_name), 'fontsize', 18, 'fontweight', 'bold')
                    end

                end                
            end    
            
            figure(1)
            set(gcf, 'resize', 'off');
            set(gcf, 'units', 'inches', 'position', [161          90        1090         667]/user_screen_resolution);
            set(gcf, 'resize', 'off');
            drawnow();
            if save_figures
                filename = ['uncertainties-nt=' num2str(num_trials, '%d') '-nr='  num2str(num_resampling_trials, '%d')];
                save_figure_to_png_svg_fig(figure_save_directory, filename, [1, 0, 0]);            
            end
            
            figure(2)
            set(gcf, 'resize', 'off');
            set(gcf, 'units', 'inches', 'position', [161          90        1090         667]/user_screen_resolution);
            set(gcf, 'resize', 'off');
            drawnow();
            if save_figures
                filename = ['snr-nt=' num2str(num_trials, '%d') '-nr='  num2str(num_resampling_trials, '%d')];
                save_figure_to_png_svg_fig(figure_save_directory, filename, [1, 0, 0]);            
            end

            return;
        end
        close all;
    end
end
% stop timer
toc

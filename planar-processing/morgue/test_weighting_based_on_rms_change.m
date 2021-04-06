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
dataset_name_array = {'PivChal03B'}; %'PivChal05B'; 'stagnation_flow'; 'Vortex_Ring'; 'Jetdata'};
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
num_trials = 1e3;

% ============================
%% resampling settings
% ============================
% resampling method names
resampling_method_names = {'remove-paired'; 'remove-random'; 'add-random'};
% resampling_method_names = {'add-random'};
num_resampling_methods = numel(resampling_method_names);
% number of particles to remove
percentage_particles_remove_array = 0.05:0.05:0.25;
num_ppr = numel(percentage_particles_remove_array);
% number of resampling trials
num_resampling_trials = 1e2; %3;
% display_resampled_images? (true/false)
display_resampled_images = 0;
% method names for plotting
resampling_method_names_plot = {'Removing Paired Particles'; 'Removing Random Particles'; 'Adding Random Particles'};
% % case name
% resampling_case_name = ['ppr' num2str(percentage_particles_remove, '%.2f') '_nrs' num2str(num_resampling_trials) '-prewindow-randpart-snr-newproc'];

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
%% weight settings
% ========================
% number of bins for weights
num_bins_w = 30;
% bins for weights
bins_w = linspace(0, 1, num_bins_w);

% ========================
%% statistics settings
% ========================
max_error_threshold = 0.2;

% ========================
%% plot settings
% ========================
% display figures?
display_figures = 1; %0;
% user screen resolution
user_screen_resolution = 113;
% save figure? (true/false)
save_figures = 0;
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
        current_write_directory = fullfile(top_write_directory, dataset_name, ['WS' num2str(window_resolution_index)], 'weights-rms-change-study');
        mkdir_c(current_write_directory);
        % file name
        filename = ['resampling-statistics-nt=' num2str(num_trials, '%d') 'nr='  num2str(num_resampling_trials, '%d') '.mat'];
        % save results to file
        load(fullfile(current_write_directory, filename));
        
        for resampling_method_index = 1:num_resampling_methods        
            % create directory to save figures
            if save_figures
                figure_save_directory = fullfile(current_write_directory, 'figures', resampling_method_names{resampling_method_index});
                mkdir_c(figure_save_directory);
            end

            wt_trials = nans(num_trials, num_individual_methods);
            h = figure;
            % ================================================
            %% loop through trials
            % ================================================
            for trial_index = 1:num_trials
                fprintf('trial_index: %d\n', trial_index);

                fprintf('resampling method index: %d\n', resampling_method_index);            
                resampling_method_name = resampling_method_names{resampling_method_index};            
                unc_ratio_rms = nans(num_individual_methods, num_ppr);
                unc_ratio_rms_rate = nans(1, num_individual_methods);

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

                        % calculate rms of the ratio
                        unc_ratio_rms(individual_method_index, particle_remove_index) = rms(unc_ratio(:, particle_remove_index), 'omitnan');

                        % name for this case
                        ppr_name_array{particle_remove_index} = num2str(round(percentage_particles_remove * 100), '%d');
                    end

                    % calculate rate of change of rms
                    fitobject = fit(percentage_particles_remove_array(1:end)', unc_ratio_rms(individual_method_index, 1:end)', 'poly1');

                    % calculate slope
                    c = coeffvalues(fitobject);
                    unc_ratio_rms_rate(individual_method_index) = c(1);

                    % ==========================
                    % plot variation of rms with particle % for just one trial
                    % ==========================
                    if trial_index == 100
                        figure(h)
                        plot(percentage_particles_remove_array*100, unc_ratio_rms(individual_method_index, :), '*', 'color', colors(individual_method_index, :))
                        hold on
                        l(individual_method_index) = plot(percentage_particles_remove_array*100, fitobject(percentage_particles_remove_array), ...
                                                            'color', [colors(individual_method_index, :), 0.5]);
                    end
                end

                % annotate figure
                if trial_index == 100
                    box off
                    xlabel('Particle Perturbation %')
                    ylabel('Uncertainty Ratio RMS')
                    legend(l, individual_method_array, 'location', 'northoutside', 'orientation', 'horizontal')
                    title(resampling_method_names_plot{resampling_method_index})

                    if save_figures
                        save_figure_to_png_svg_fig(figure_save_directory, 'rms-rate', [1, 0, 0])
                    end
                end

                % calculate weights
                wt = 1 ./ abs(unc_ratio_rms_rate);
                wt_trials(trial_index, :) = wt ./ sum(wt);
            end
            
            % ==========================
            % calculate and plot pdf of weights
            % ==========================
            pdf_w = nans(num_individual_methods, num_bins_w-1);
            figure
            for individual_method_index = 1:num_individual_methods
                pdf_w(individual_method_index, :) = histcounts(wt_trials(:, individual_method_index), bins_w, 'normalization', 'pdf');
                
                % subplot(num_individual_methods, 1, individual_method_index)
                plot(bins_w(1:end-1), pdf_w(individual_method_index, :), 'color', colors(1, :), 'linestyle', line_symbols{individual_method_index})                        
                hold on
                % title(individual_method_array{individual_method_index});

                % if individual_method_index < num_individual_methods
                %     set(gca, 'xticklabel', [])
                % else
                %     xlabel('Weights')
                % end
            end
            box off
            xlim([0 1]);
            xlabel('Weights')
            ylabel('Probability Density')
            legend(individual_method_array, 'Location', 'northoutside', 'Orientation', 'horizontal')
            title(resampling_method_names_plot{resampling_method_index})

            fig = gcf;
            set(gcf, 'resize', 'off');
            set(gcf, 'units', 'inches', 'position', [440   262   494   512]/user_screen_resolution);
            set(gcf, 'resize', 'off');

            if save_figures
                save_figure_to_png_svg_fig(figure_save_directory, 'pdf-weights', [1, 0, 0]);
            end

            % ==========================
            % extract errors and uncertainties
            % ==========================
            err_U_trials = nans(num_trials, 1);        
            errors = errors_all{window_resolution_index, dataset_index}.err_U;
            snapshot_index = 1;
            for trial_index = 1:num_trials
                err_U_trials(trial_index) = errors(r_trials(trial_index), c_trials(trial_index), snapshot_index);            
            end

            err_rms = rms(err_U_trials, 'omitnan');
            % err_rms = std(err_U_trials, [], 'omitnan');
            % calculate confidence interval for the error
            confint_err = prctile(err_U_trials, 84) - prctile(err_U_trials, 16);
            confint_err_abs = prctile(abs(err_U_trials), 84) - prctile(abs(err_U_trials), 16);

            % ==========================
            % extract individual uncertainty
            % ==========================
            unc_indiv = nans(num_trials, num_individual_methods);
            results = results_all{window_resolution_index, dataset_index}{snapshot_index};
            for trial_index = 1:num_trials
                unc_current = extract_planar_uncertainties(results.uncertainty2D, r_trials(trial_index), c_trials(trial_index));
                for individual_method_index = 1:num_individual_methods
                    method_name = [lower(individual_method_array{individual_method_index}) 'x'];
                    unc_indiv(trial_index, individual_method_index) = unc_current.(method_name);                
                end
            end

            sigma_rms_indv = rms(unc_indiv, 1);

            % ==========================
            % calculate combined uncertainty
            % ==========================
            unc_comb = sum(unc_indiv  .* wt_trials, 2); 
            sigma_rms_comb = rms(unc_comb, 'omitnan');

            % ==========================
            % plot pdfs of errors and uncertainties
            % ==========================
            figure
            violins = violinplot([unc_indiv, abs(err_U_trials), unc_comb], {individual_method_array{:}, 'Error', 'RMS-Rate'} , ...
                'showdata', false, 'shownotches', true, 'edgecolor', [1, 1, 1], 'ViolinAlpha', 0.5);

            % face colors
            color_all = cell(1, numel(violins));

            % add lines corresponding to rms
            for violin_index = 1:numel(violins)
                % x co-ordinates of current violin plot
                x = violins(violin_index).ViolinPlot.XData;

                % individual uncertainty methods
                if violin_index <= num_individual_methods
                    violins(violin_index).ViolinColor = colors(1, :);
                    % violins(violin_index).ViolinColor = colors_blue{violin_index-1};
                    y = sigma_rms_indv(violin_index);
                % error
                elseif violin_index == num_individual_methods+1
                    violins(violin_index).ViolinColor = [0, 0, 0];
                    y = err_rms;
                    % x = [violins(violin_index).ViolinPlot.XData; violins(numel(violins)).ViolinPlot.XData];
                    plot([min(violins(1).ViolinPlot.XData), max(violins(numel(violins)).ViolinPlot.XData)], y*[1, 1], '--', 'color', [0, 0, 0, 0.2])

                    % plot confidence intervals for the error
                    plot([min(violins(1).ViolinPlot.XData), max(violins(numel(violins)).ViolinPlot.XData)], confint_err*[1, 1], '-*', 'color', [0, 0, 0, 0.2])
                    plot([min(violins(1).ViolinPlot.XData), max(violins(numel(violins)).ViolinPlot.XData)], confint_err_abs*[1, 1], '-o', 'color', [0, 0, 0, 0.2])
    
                % combined uncertainty methods
                elseif violin_index > num_individual_methods + 1
                    violins(violin_index).ViolinColor = colors(2, :);
                    % violins(violin_index).ViolinColor = colors_red{violin_index-4};
                    y = sigma_rms_comb(violin_index - (num_individual_methods + 1));
                end

                %% adjust violin properties
                color_all{violin_index} = violins(violin_index).ViolinColor;
                violins(violin_index).ViolinAlpha = 0.25;
                violins(violin_index).BoxColor = violins(violin_index).ViolinColor;
                violins(violin_index).BoxWidth = 0.005;

                %% plot rms value
                plot([min(x), max(x)], y*[1, 1], 'color', [violins(violin_index).ViolinColor, violins(violin_index).ViolinAlpha+0.2], 'linewidth', 3)
            end

            % annotate y axis
            ylabel('Error/Uncertainty (pix.)', 'fontsize', 16)
            % adjust limits
            ylim([0 max_error_threshold])
            % ylim([0 0.15])

            pause(0.1);
            % turn off x axis line
            ax = gca;
            ax.XAxis.Axle.Visible = 'off';
            ax.XAxis.TickLength = [0 0];
            % turn off y axis line
            ax.YAxis.Axle.Visible = 'off';
            % ax.YAxis.TickLength = [0, 0];
            ax.YAxis.TickLength = [0.005 0.005];

            title(resampling_method_names_plot{resampling_method_index})            
            % adjust figure position
            set(gcf, 'resize', 'off');
            drawnow();
            set(gcf, 'units', 'inches', 'Position', [352   526   895   392]/user_screen_resolution)
            drawnow();

            if save_figures
                save_figure_to_png_svg_fig(figure_save_directory, 'violin', [1, 0, 0]);
            end

            return;
        end
        return;
    end
end
% stop timer
toc

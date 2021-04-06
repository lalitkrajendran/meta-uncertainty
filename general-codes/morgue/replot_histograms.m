clear
close all
clc

addpath('/scratch/shannon/c/aether/Projects/BOS/error-analysis/analysis/src/piv-image-generation');
addpath('/scratch/shannon/c/aether/Projects/BOS/error-analysis/analysis/src/piv-image-generation/jobfiles');
addpath prana-uncertainty-average-dc-new-im-2/
addpath(genpath('/scratch/shannon/c/aether/Projects/BOS/general-codes/matlab-codes/'));
setup_default_settings;
% addpath('/scratch/shannon/c/aether/Projects/BOS/uncertainty/analysis/src/prana-uncertainty-average-dc-new-im-2/');

logical_string = {'False'; 'True'};

%% read/write settings

% top level directory for saving images
top_image_directory = fullfile('/scratch/shannon/a/aether/Projects/PIV_PTV_Uncertainty_Quantification/uncertainty-combination/analysis/data/images/local-averaging/');
if ~exist(top_image_directory, 'dir')
    mkdir(top_image_directory);
end

% top results directory
top_results_directory = fullfile('/scratch/shannon/a/aether/Projects/PIV_PTV_Uncertainty_Quantification/uncertainty-combination/analysis/results/local-averaging/');
if ~exist(top_results_directory, 'dir')
    mkdir(top_results_directory);
end
% save tiff images? (true/false)
save_images_to_tiff = false;

%% image properties
% mean displacement (pix.)
delta_x = 0;
% standard deviation of displacement (pix.)
displacement_noise_std = 0.1;
% number of image pairs to generate
num_image_pairs = 1e3;
% height of image (pix.)
image_height_array = [32, 64];
% % width of image (pix.)
% image_width = 32;
% particle diameter (pix.)
d_p = 3;
% array of image noise levels
image_noise_std_array = [0, 5]; 
% % particle overlap? (true/false)
% particle_overlap = true;

%% processing settings

% uncertainty flags for the various methods (0 = False, 1 = True)
uncertainty.ppruncertainty = 0;
uncertainty.miuncertainty = 0;
uncertainty.mcuncertainty = 1;

% minimum uncertainty level
min_uncertainty_threshold = 1e-4;
max_uncertainty_threshold = 0.1;
%% bootstrapping settings

% number of particles to remove
% num_particles_remove_array = [5, 10, 15];
percentage_particles_remove = 0.2;
% number of bootstrap trials
num_bootstrap_trials = 1e3;

%% plot settings

% bins for histogram
% bins = linspace(0, displacement_noise_std*2, round(sqrt(num_image_pairs)));
bins = linspace(min_uncertainty_threshold, max_uncertainty_threshold, round(sqrt(num_image_pairs)));

%% analysis
for image_height = image_height_array
    fprintf('image_height: %d\n', image_height);
    image_width = image_height;
    for particle_overlap = [false, true]
        fprintf('particle overlap: %s\n', logical_string{particle_overlap + 1});
        for image_noise_std = image_noise_std_array
            fprintf('noise: %d\n', image_noise_std);
            %% create directories
            % directory where images for this case will be saved
            % current_image_directory = fullfile(top_image_directory, ['delta_x=' num2str(delta_x, '%.1f') 'pix_std=' num2str(displacement_noise_std, '%.1f') 'pix'], [num2str(image_height) 'x' num2str(image_width)], ['noise' num2str(image_noise_std, '%02d')]);
            current_image_directory = fullfile(top_image_directory, ['overlap=' logical_string{particle_overlap + 1}], ['delta_x=' num2str(delta_x, '%.1f') 'pix_std=' num2str(displacement_noise_std, '%.1f') 'pix_imnoise=' num2str(image_noise_std, '%02d')], [num2str(image_height) 'x' num2str(image_width)]);
            if ~exist(current_image_directory, 'dir')
                mkdir(current_image_directory);
            end

            %% estimate disparity from all image pairs
            fprintf('Estimating disparity by monte-carlo (global)\n');
            % set grid point locations
            X = image_width/2;
            Y = image_height/2;

            % set displacements
            U = delta_x;
            V = delta_x;

            % load particle positions
            particle_positions = load(fullfile(current_image_directory, 'parameters', ['particle_positions_lin_h' num2str(image_height)   '_w' num2str(image_width) '_seg_000001_' num2str(num_image_pairs, '%06d') '.mat']));

            % calculate number of particles in the interrogation window
            num_particles = numel(particle_positions.ParticlePositions(1).X1);

            %% estimate disparity histogram for all image pairs by bootstrapping

            num_particles_remove = round(percentage_particles_remove * num_particles);
            fprintf('num particles remove: %d\n', num_particles_remove);

            % directory where results for this case will be saved
            current_results_directory = fullfile(top_results_directory, ['overlap=' logical_string{particle_overlap + 1}], ['delta_x=' num2str(delta_x, '%.1f') 'pix_std=' num2str(displacement_noise_std, '%.1f') 'pix_imnoise=' num2str(image_noise_std, '%02d')], [num2str(image_height) 'x' num2str(image_width)], ['ntot' num2str(num_particles, '%02d') '_npr' num2str(num_particles_remove, '%02d') '_nbs' num2str(num_bootstrap_trials)]);
            if ~exist(current_results_directory, 'dir')
                mkdir(current_results_directory);
            end

            fprintf('Estimating disparity by bootstrapping (local)\n');

            %% load results

            % load uncertainties from individual schemes
            load(fullfile(current_results_directory, 'sigma_individual.mat'));

            % load weighted uncertainties
            load(fullfile(current_results_directory, 'sigma_combined.mat'));

            % load weights
            load(fullfile(current_results_directory, 'weights.mat'));

            % load pdfs
            load(fullfile(current_results_directory, 'histograms.mat'));

            %% plot histogram of uncertainties

            % directory to save figures
            figure_save_directory = fullfile(current_results_directory, 'figures');
            mkdir_c(figure_save_directory);
            
            y_max = 150;

            figure

            l_im = plot(bins(1:end-1), N_im_global, 'c');
            hold on
            plot(sigma_rms_im_global * ones(1,100), linspace(0, y_max, 100), 'c')

            l_mc = plot(bins(1:end-1), N_mc_global, 'r');
            plot(sigma_rms_mc_global * ones(1,100), linspace(0, y_max, 100), 'r')

            l_cs = plot(bins(1:end-1), N_cs_global, 'm');
            plot(sigma_rms_cs_global * ones(1,100), linspace(0, y_max, 100), 'm')

            l_wt_none = plot(bins(1:end-1), N_weighted_none, 'g');
            plot(sigma_rms_weighted_none * ones(1,100), linspace(0, y_max, 100), 'g')

            l_wt_global = plot(bins(1:end-1), N_weighted_global, 'g--');
            plot(sigma_rms_weighted_global * ones(1,100), linspace(0, y_max, 100), 'g--')

            l_wt_local = plot(bins(1:end-1), N_weighted_local, 'g:');
            plot(sigma_rms_weighted_local * ones(1,100), linspace(0, y_max, 100), 'g:')

            l_true = plot(sigma_rms_true * ones(1,100), linspace(0, y_max, 100), 'k');

            xlim([0 0.1])
            ylim([0 y_max])
            legend([l_im, l_mc, l_cs, l_wt_none, l_wt_global, l_wt_local, l_true], ...
                {'\sigma_{IM}', '\sigma_{MC}', '\sigma_{CS}', 'Unweighted', 'Global', 'Local', 'RMS, True'}, 'fontsize', 12, ...
                'location', 'eastoutside')
            
            xlabel('Uncertainty (pix.)')
            ylabel('Probability Density', 'fontsize', 12)

            set(gcf, 'Position', [360   573   628   327])
            save_figure_to_png_eps_fig(figure_save_directory, 'histograms-comparison', [1, 0, 0]);
            %%
            close all
        end
    end
end
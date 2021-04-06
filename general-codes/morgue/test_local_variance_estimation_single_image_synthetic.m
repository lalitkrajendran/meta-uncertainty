clear
close all
clc

addpath('/scratch/shannon/c/aether/Projects/BOS/error-analysis/analysis/src/piv-image-generation');
addpath('/scratch/shannon/c/aether/Projects/BOS/error-analysis/analysis/src/piv-image-generation/jobfiles');
addpath prana-uncertainty-average-dc-new-im-2/
addpath(genpath('/scratch/shannon/c/aether/Projects/BOS/general-codes/matlab-codes/'));

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
num_image_pairs = 2; %1e4;
% height of image (pix.)
image_height = 64;
% width of image (pix.)
image_width = 64;
% particle diameter (pix.)
d_p = 3;
% array of image noise levels
image_noise_std_array = [0, 5, 10];
% particle overlap? (true/false)
particle_overlap = true;

%% processing settings

% uncertainty flags for the various methods (0 = False, 1 = True)
uncertainty.ppruncertainty = 0;
uncertainty.miuncertainty = 0;
uncertainty.mcuncertainty = 1;

%% bootstrapping settings

% number of particles to remove
num_particles_remove = 10;
% number of bootstrap trials
num_bootstrap_trials = 1e4;

%% plot settings

% bins for histogram
bins = linspace(0, displacement_noise_std*2, round(sqrt(num_image_pairs)));

%% analysis
for image_noise_std = image_noise_std_array
    fprintf('noise: %d\n', image_noise_std);
    %% create directories
    % directory where images for this case will be saved
    % current_image_directory = fullfile(top_image_directory, ['delta_x=' num2str(delta_x, '%.1f') 'pix_std=' num2str(displacement_noise_std, '%.1f') 'pix'], [num2str(image_height) 'x' num2str(image_width)], ['noise' num2str(image_noise_std, '%02d')]);
    current_image_directory = fullfile(top_image_directory, ['overlap=' logical_string{particle_overlap + 1}], ['delta_x=' num2str(delta_x, '%.1f') 'pix_std=' num2str(displacement_noise_std, '%.1f') 'pix_imnoise=' num2str(image_noise_std, '%02d')], [num2str(image_height) 'x' num2str(image_width)]);
    if ~exist(current_image_directory, 'dir')
        mkdir(current_image_directory);
    end

    % directory where results for this case will be saved
    current_results_directory = fullfile(top_results_directory, ['overlap=' logical_string{particle_overlap + 1}], ['delta_x=' num2str(delta_x, '%.1f') 'pix_std=' num2str(displacement_noise_std, '%.1f') 'pix_imnoise=' num2str(image_noise_std, '%02d')], [num2str(image_height) 'x' num2str(image_width)], ['npr' num2str(num_particles_remove) '_nbs' num2str(num_bootstrap_trials)]);
    if ~exist(current_results_directory, 'dir')
        mkdir(current_results_directory);
    end

    %% generate image matrices

    % generate job file
    % job_file = MonteCarloImageGenerationJobFile_displacement_noise(displacement_noise_std, delta_x, num_image_pairs, image_height, image_width, d_p, top_image_directory);
    job_file = MonteCarloImageGenerationJobFile_general(displacement_noise_std, delta_x, num_image_pairs, image_height, image_width, d_p, image_noise_std*0.01, fullfile(top_image_directory, ['overlap=' logical_string{particle_overlap + 1}]));
    % generate images
    [image_matrix_1, image_matrix_2] = generateMonteCarloImageSet(job_file, particle_overlap);

    %% convert matrix to images and save
    if save_images_to_tiff
        % directory to save images
        image_write_directory = fullfile(current_image_directory, 'image');
        if ~exist(image_write_directory, 'dir')
            mkdir(image_write_directory);
        end

        % loop through all image pairs
        for image_pair_index = 1:num_image_pairs
            fprintf('image_pair_index: %d\n', image_pair_index);

            % filename for first image of the pair
            image_filename = ['img_' num2str(2*(image_pair_index-1) + 1, '%04d') '.tif'];
            % save im1
            imwrite(image_matrix_1(:,:,image_pair_index), fullfile(image_write_directory, image_filename));

            % filename for im2
            image_filename = ['img_' num2str(2*image_pair_index, '%04d') '.tif'];
            % save im2
            imwrite(image_matrix_2(:,:,image_pair_index), fullfile(image_write_directory, image_filename));
        end
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

    % initialize arrays
    sigma_x_im_global = nans(1, num_image_pairs);
    sigma_y_im_global = nans(1, num_image_pairs);
    sigma_x_mc_global = nans(1, num_image_pairs);
    sigma_y_mc_global = nans(1, num_image_pairs);
    sigma_x_cs_global = nans(1, num_image_pairs);
    sigma_y_cs_global = nans(1, num_image_pairs);

    % loop through image pairs and calculate disparity
    for image_pair_index = 1:num_image_pairs
        % extract image pair
        im1 = double(image_matrix_1(:, :, image_pair_index));
        im2 = double(image_matrix_2(:, :, image_pair_index));

        %% calculate true disparity

        % extract particle positions for the first image
        x1 = particle_positions.ParticlePositions(image_pair_index).X1;
        y1 = particle_positions.ParticlePositions(image_pair_index).Y1;

        % extract particle positions for the second image
        x2 = particle_positions.ParticlePositions(image_pair_index).X2;
        y2 = particle_positions.ParticlePositions(image_pair_index).Y2;    

        % calculate disparity along x
        disparity_x_true = std(x2 - x1) / sqrt(numel(x1));
        disparity_y_true = std(y2 - y1) / sqrt(numel(y1));

        %% calculate uncertainty using IM
        [sigma_x_im_global(image_pair_index), sigma_y_im_global(image_pair_index), Npartu, Npartv] = original_particle_disparity_no_deform(im1, im2, X, Y, U, V, image_height);

        %% calculate uncertainty using MC

        [~,~,~,~,~,~,~,uncertainty2D,~] = PIVwindowed(im1,im2,'SCC',[image_height, image_width], [image_height, image_width; image_height, image_width], 0, [2.8, 2.8], 0, 3, 0, 0, 0, X, Y, uncertainty, 1, 0, 0);

        sigma_x_mc_global(image_pair_index) = sqrt(uncertainty2D.biasx.^2+(uncertainty2D.Ixx.^2)./uncertainty2D.Neff);
        sigma_y_mc_global(image_pair_index) = sqrt(uncertainty2D.biasy.^2+(uncertainty2D.Iyy.^2)./uncertainty2D.Neff);

        %% calculate uncertainty using CS
        [sigma_x_cs_global(image_pair_index), sigma_y_cs_global(image_pair_index)] = correlation_statistics(im1,im2,[image_height, image_width],[image_height, image_width], X, Y);

    end

    %% calcualte standard deviation of disparity from monte carlo

    std_sigma_x_im_global = std(sigma_x_im_global, 'omitnan');
    std_sigma_y_im_global = std(sigma_y_im_global, 'omitnan');
    std_sigma_x_mc_global = std(sigma_x_mc_global, 'omitnan');
    std_sigma_y_mc_global = std(sigma_y_mc_global, 'omitnan');
    std_sigma_x_cs_global = std(sigma_x_cs_global, 'omitnan');
    std_sigma_y_cs_global = std(sigma_y_cs_global, 'omitnan');

    %% calculate global weights
    w_x_im_global = 1/std_sigma_x_im_global^2/(1/std_sigma_x_im_global^2 + 1/std_sigma_x_mc_global^2 + 1/std_sigma_x_cs_global^2);
    w_y_im_global = 1/std_sigma_y_im_global^2/(1/std_sigma_y_im_global^2 + 1/std_sigma_y_mc_global^2 + 1/std_sigma_y_cs_global^2);

    w_x_mc_global = 1/std_sigma_x_mc_global^2/(1/std_sigma_x_im_global^2 + 1/std_sigma_x_mc_global^2 + 1/std_sigma_x_cs_global^2);
    w_y_mc_global = 1/std_sigma_y_mc_global^2/(1/std_sigma_y_im_global^2 + 1/std_sigma_y_mc_global^2 + 1/std_sigma_y_cs_global^2);

    w_x_cs_global = 1 - (w_x_im_global + w_x_mc_global);
    w_y_cs_global = 1 - (w_y_im_global + w_y_mc_global);

    %% estimate disparity histogram for a single image pair by bootstrapping
    fprintf('Estimating disparity by bootstrapping (local)\n');

    % set image pair index
    image_pair_index = 1;

    % extract current image
    im1 = double(image_matrix_1(:, :, image_pair_index));
    im2 = double(image_matrix_2(:, :, image_pair_index));

    % extract particle positions for the first image
    x1 = particle_positions.ParticlePositions(image_pair_index).X1;
    y1 = particle_positions.ParticlePositions(image_pair_index).Y1;

    % extract particle positions for the second image
    x2 = particle_positions.ParticlePositions(image_pair_index).X2;
    y2 = particle_positions.ParticlePositions(image_pair_index).Y2;

    % initialize arrays
    sigma_x_im_local = nans(1, num_bootstrap_trials);
    sigma_y_im_local = nans(1, num_bootstrap_trials);
    sigma_x_mc_local = nans(1, num_bootstrap_trials);
    sigma_y_mc_local = nans(1, num_bootstrap_trials);
    sigma_x_cs_local = nans(1, num_bootstrap_trials);
    sigma_y_cs_local = nans(1, num_bootstrap_trials);

    parfor trial_index = 1:num_bootstrap_trials

        % randomly sample particles for removal
        particle_sample = randsample(numel(x1), num_particles_remove);

        %% remove the pixels surrounding these particles from the images
        im1_temp = im1;
        im2_temp = im2;

        for sample_index = 1:num_particles_remove

            % index of the particle that is to be removed
            particle_index = particle_sample(sample_index);

            % x, y location of the particle to be removed
            x1_current = x1(particle_index);
            y1_current = y1(particle_index);
            x2_current = x2(particle_index);
            y2_current = y2(particle_index);

            % integer locations on the image
            r1_array = (-d_p:d_p) + image_height - round(y1_current) + 1;
            c1_array = (-d_p:d_p) + round(x1_current);
            r2_array = (-d_p:d_p) + image_height - round(y2_current) + 1;
            c2_array = (-d_p:d_p) + round(x2_current);

            % loop through locations and set pixel intensities to zero
            for r_index = 1:numel(r1_array)
                for c_index = 1:numel(c1_array)

                    % row, col locations for the current particle in im1
                    r1_current = r1_array(r_index);
                    c1_current = c1_array(c_index);

                    % zero out intensity map
                    if check_valid_index(image_height, image_width, r1_current, c1_current)
                        im1_temp(r1_current, c1_current) = 0;
                    end

                    % row, col locations for the current particle in im2
                    r2_current = r2_array(r_index);
                    c2_current = c2_array(c_index);

                    % zero out intensity map
                    if check_valid_index(image_height, image_width, r2_current, c2_current)
                        im2_temp(r2_current, c2_current) = 0;
                    end                
                end
            end
        end

        %% calculate uncertainty using IM
        [sigma_x_im_local(trial_index), sigma_y_im_local(trial_index), ~, ~] = original_particle_disparity_no_deform(im1_temp, im2_temp, X, Y, U, V, image_height);

        %% calculate uncertainty using MC
        [~,~,~,~,~,~,~,uncertainty2D,~] = PIVwindowed(im1_temp,im2_temp,'SCC',[image_height, image_width], [image_height, image_width; image_height, image_width], 0, [2.8, 2.8], 0, 3, 0, 0, 0, X, Y, uncertainty, 1, 0, 0);

        sigma_x_mc_local(trial_index) = sqrt(uncertainty2D.biasx.^2+(uncertainty2D.Ixx.^2)./uncertainty2D.Neff);
        sigma_y_mc_local(trial_index) = sqrt(uncertainty2D.biasy.^2+(uncertainty2D.Iyy.^2)./uncertainty2D.Neff);

        %% calculate uncertainty using CS
        [sigma_x_cs_local(trial_index), sigma_y_cs_local(trial_index)] = correlation_statistics(im1_temp,im2_temp,[image_height, image_width],[image_height, image_width], X, Y);

    end

    %% calculate standard deviation of disparity from bootstrapping
    std_sigma_x_im_local = std(sigma_x_im_local, 'omitnan');
    std_sigma_y_im_local = std(sigma_y_im_local, 'omitnan');
    std_sigma_x_mc_local = std(sigma_x_mc_local, 'omitnan');
    std_sigma_y_mc_local = std(sigma_y_mc_local, 'omitnan');
    std_sigma_x_cs_local = std(sigma_x_cs_local, 'omitnan');
    std_sigma_y_cs_local = std(sigma_y_cs_local, 'omitnan');

    %% calculate local weights

    w_x_im_local = 1/std_sigma_x_im_local^2/(1/std_sigma_x_im_local^2 + 1/std_sigma_x_mc_local^2 + 1/std_sigma_x_cs_local^2);
    w_y_im_local = 1/std_sigma_y_im_local^2/(1/std_sigma_y_im_local^2 + 1/std_sigma_y_mc_local^2 + 1/std_sigma_y_cs_local^2);

    w_x_mc_local = 1/std_sigma_x_mc_local^2/(1/std_sigma_x_im_local^2 + 1/std_sigma_x_mc_local^2 + 1/std_sigma_x_cs_local^2);
    w_y_mc_local = 1/std_sigma_y_mc_local^2/(1/std_sigma_y_im_local^2 + 1/std_sigma_y_mc_local^2 + 1/std_sigma_y_cs_local^2);

    w_x_cs_local = 1 - (w_x_im_local + w_x_mc_local);
    w_y_cs_local = 1 - (w_y_im_local + w_y_mc_local);

    %% compare standard deviations of the disparity estimaets    
    fprintf('writing statistics to file\n');
    
    fileID = fopen(fullfile(current_results_directory, 'statistics.txt'), 'w');
    
    fprintf(fileID, '-------------------------------\n');
    fprintf(fileID, 'Image Matching\n');
    fprintf(fileID, '-------------------------------\n');
    fprintf(fileID, 'Std disparity, Global: %.3f, %.3f pix.\n', std_sigma_x_im_global, std_sigma_x_im_global);
    fprintf(fileID, 'Std disparity, BS: %.3f, %.3f pix.\n', std_sigma_x_im_local, std_sigma_y_im_local);
    fprintf(fileID, 'Ratio: %.3f, %.3f pix.\n', std_sigma_x_im_local/std_sigma_x_im_global, std_sigma_y_im_local/std_sigma_y_im_global);
    
    fprintf(fileID, '-------------------------------\n');
    fprintf(fileID, 'Moment of Correlation\n');
    fprintf(fileID, '-------------------------------\n');
    fprintf(fileID, 'Std disparity, Global: %.3f, %.3f pix.\n', std_sigma_x_mc_global, std_sigma_x_mc_global);
    fprintf(fileID, 'Std disparity, BS: %.3f, %.3f pix.\n', std_sigma_x_mc_local, std_sigma_y_mc_local);
    fprintf(fileID, 'Ratio: %.3f, %.3f pix.\n', std_sigma_x_mc_local/std_sigma_x_mc_global, std_sigma_y_mc_local/std_sigma_y_mc_global);

    fprintf(fileID, '-------------------------------\n');
    fprintf(fileID, 'Correlation Statistics\n');
    fprintf(fileID, '-------------------------------\n');
    fprintf(fileID, 'Std disparity, Global: %.3f, %.3f pix.\n', std_sigma_x_cs_global, std_sigma_x_cs_global);
    fprintf(fileID, 'Std disparity, BS: %.3f, %.3f pix.\n', std_sigma_x_cs_local, std_sigma_y_cs_local);
    fprintf(fileID, 'Ratio: %.3f, %.3f pix.\n', std_sigma_x_cs_local/std_sigma_x_cs_global, std_sigma_y_cs_local/std_sigma_y_cs_global);

    fclose(fileID);

    %% compare weights
    fprintf('writing weights to file\n');
    
    fileID = fopen(fullfile(current_results_directory, 'weights.txt'), 'w');

    fprintf(fileID, '-------------------------------\n');
    fprintf(fileID, 'Weights\n');
    fprintf(fileID, '-------------------------------\n');
    fprintf(fileID, 'Global-x: IM: %.3f, \tMC:%.3f, CS:%.3f \n', w_x_im_global, w_x_mc_global, w_x_cs_global);
    fprintf(fileID, 'Local-x: IM: %.3f, MC:%.3f, CS:%.3f \n', w_x_im_local, w_x_mc_local, w_x_cs_local);
    fprintf(fileID, 'Global-y: IM: %.3f, MC:%.3f, CS:%.3f \n', w_y_im_global, w_y_mc_global, w_y_cs_global);
    fprintf(fileID, 'Local-y: IM: %.3f, MC:%.3f, CS:%.3f \n', w_y_im_local, w_y_mc_local, w_y_cs_local);

    fclose(fileID);
    
    %% plot disparity histograms

    % directory to save figures
    figure_save_directory = fullfile(current_results_directory, 'figures');
    mkdir_c(figure_save_directory);

    y_max = 200;

    % -------------------------------
    % Image Matching
    % -------------------------------

    figure
    subplot(2, 2, 1)
    histogram(sigma_x_im_global, bins, 'normalization', 'pdf')
    xlim([0, displacement_noise_std])
    ylim([0 y_max])
    xlabel('pix.')
    title({'DX, Global';['\sigma = ' num2str(std_sigma_x_im_global, '%.3f') ' pix.']})

    subplot(2, 2, 2)
    histogram(sigma_y_im_global, bins, 'normalization', 'pdf')
    xlim([0, displacement_noise_std])
    ylim([0 y_max])
    xlabel('pix.')
    title({'DY, Global';['\sigma = ' num2str(std_sigma_y_im_global, '%.3f') ' pix.']})

    subplot(2, 2, 3)
    histogram(sigma_x_im_local, bins, 'normalization', 'pdf')
    xlim([0, displacement_noise_std])
    ylim([0 y_max])
    xlabel('pix.')
    title({'DX, BS';['\sigma = ' num2str(std_sigma_x_im_local, '%.3f') ' pix.']})

    subplot(2, 2, 4)
    histogram(sigma_y_im_local, bins, 'normalization', 'pdf')
    xlim([0, displacement_noise_std])
    ylim([0 y_max])
    xlabel('pix.')
    title({'DY, BS';['\sigma = ' num2str(std_sigma_y_im_local, '%.3f') ' pix.']})

    set(gcf, 'Position', [360   362   697   538])

    save_figure_to_png_eps_fig(figure_save_directory, 'histograms-im', [1, 0, 0]); 

    % -------------------------------
    % Moment of Correlation
    % -------------------------------

    figure
    subplot(2, 2, 1)
    histogram(sigma_x_mc_global, bins, 'normalization', 'pdf')
    xlim([0, displacement_noise_std])
    ylim([0 y_max])
    xlabel('pix.')
    title({'DX, Global';['\sigma = ' num2str(std_sigma_x_mc_global, '%.3f') ' pix.']})

    subplot(2, 2, 2)
    histogram(sigma_y_mc_global, bins, 'normalization', 'pdf')
    xlim([0, displacement_noise_std])
    ylim([0 y_max])
    xlabel('pix.')
    title({'DY, Global';['\sigma = ' num2str(std_sigma_y_mc_global, '%.3f') ' pix.']})

    subplot(2, 2, 3)
    histogram(sigma_x_mc_local, bins, 'normalization', 'pdf')
    xlim([0, displacement_noise_std])
    ylim([0 y_max])
    xlabel('pix.')
    title({'DX, BS';['\sigma = ' num2str(std_sigma_x_mc_local, '%.3f') ' pix.']})

    subplot(2, 2, 4)
    histogram(sigma_y_mc_local, bins, 'normalization', 'pdf')
    xlim([0, displacement_noise_std])
    ylim([0 y_max])
    xlabel('pix.')
    title({'DY, BS';['\sigma = ' num2str(std_sigma_y_mc_local, '%.3f') ' pix.']})

    set(gcf, 'Position', [360   362   697   538])

    save_figure_to_png_eps_fig(figure_save_directory, 'histograms-mc', [1, 0, 0]); 

    %%
    % -------------------------------
    % Correlation Statistics
    % -------------------------------
    figure
    subplot(2, 2, 1)
    histogram(sigma_x_cs_global, bins, 'normalization', 'pdf')
    xlim([0, displacement_noise_std])
    ylim([0 y_max])
    xlabel('pix.')
    title({'DX, Global';['\sigma = ' num2str(std_sigma_x_cs_global, '%.3f') ' pix.']})

    subplot(2, 2, 2)
    histogram(sigma_y_cs_global, bins, 'normalization', 'pdf')
    xlim([0, displacement_noise_std])
    ylim([0 y_max])
    xlabel('pix.')
    title({'DY, Global';['\sigma = ' num2str(std_sigma_y_cs_global, '%.3f') ' pix.']})

    subplot(2, 2, 3)
    histogram(sigma_x_cs_local, bins, 'normalization', 'pdf')
    xlim([0, displacement_noise_std])
    ylim([0 y_max])
    xlabel('pix.')
    title({'DX, BS';['\sigma = ' num2str(std_sigma_x_cs_local, '%.3f') ' pix.']})

    subplot(2, 2, 4)
    histogram(sigma_y_cs_local, bins, 'normalization', 'pdf')
    xlim([0, displacement_noise_std])
    ylim([0 y_max])
    xlabel('pix.')
    title({'DY, BS';['\sigma = ' num2str(std_sigma_y_cs_local, '%.3f') ' pix.']})

    set(gcf, 'Position', [360   362   697   538])

    save_figure_to_png_eps_fig(figure_save_directory, 'histograms-cs', [1, 0, 0]);
    %%
    close all
end
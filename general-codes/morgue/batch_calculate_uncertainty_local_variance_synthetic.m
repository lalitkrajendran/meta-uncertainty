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
num_image_pairs = 1; %e3;
% height of image (pix.)
image_height_array = [32, 64];
% % width of image (pix.)
% image_width = 32;
% particle diameter (pix.)
d_p = 3;
% array of image noise levels
image_noise_std_array = [5]; 
% % particle overlap? (true/false)
% particle_overlap = true;

%% processing settings

% uncertainty flags for the various methods (0 = False, 1 = True)
uncertainty_flags.ppruncertainty = 0;
uncertainty_flags.miuncertainty = 0;
uncertainty_flags.mcuncertainty = 1;

% minimum uncertainty level
min_uncertainty_threshold = 1e-4;
max_uncertainty_threshold = 0.1;

%% bootstrapping settings

% method to obtain particle position ('true', 'id')
particle_position_calc_method = 'id';
% intensity threshold for particle identification
intensity_threshold = 1e3;

% particle sizing settings
sizeprops = struct;
sizeprops.method = 'IWC';
sizeprops.p_area = 0.5 * d_p^2;
sizeprops.sigma = 4;
sizeprops.errors = 0;
sizeprops.save_dir = [];

% number of particles to remove
% num_particles_remove_array = [5, 10, 15];
percentage_particles_remove = 0.2;
% number of bootstrap trials
num_bootstrap_trials = 10; %e3;

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

            % calculate number of particles in the interrogation window
            num_particles = numel(particle_positions.ParticlePositions(1).X1);

            % initialize arrays
            sigma_x_im_global = nans(1, num_image_pairs);
            sigma_y_im_global = nans(1, num_image_pairs);
            sigma_x_mc_global = nans(1, num_image_pairs);
            sigma_y_mc_global = nans(1, num_image_pairs);
            sigma_x_cs_global = nans(1, num_image_pairs);
            sigma_y_cs_global = nans(1, num_image_pairs);

            % loop through image pairs and calculate disparity
            parfor image_pair_index = 1:num_image_pairs
    %             fprintf('image pair: %d\n', image_pair_index);
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

                [~,~,~,~,~,~,~,uncertainty2D,~] = PIVwindowed(im1,im2,'SCC',[image_height, image_width], [image_height, image_width; image_height, image_width], 0, [2.8, 2.8], 0, 3, 0, 0, 0, X, Y, uncertainty_flags, 1, 0, 0);

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

            %% calculate unweighted uncertainty

            % calculate unweighted uncertainty
            sigma_x_weighted_none = 1/3 * sigma_x_im_global + 1/3 * sigma_x_mc_global + 1/3 * sigma_x_cs_global;
            sigma_y_weighted_none = 1/3 * sigma_y_im_global + 1/3 * sigma_y_mc_global + 1/3 * sigma_y_cs_global;        

            % calculate globally weighted uncertainty
            sigma_x_weighted_global = w_x_im_global * sigma_x_im_global + ...
                w_x_mc_global * sigma_x_mc_global + ...
                w_x_cs_global * sigma_x_cs_global;
            sigma_y_weighted_global = w_y_im_global * sigma_y_im_global + ...
                w_y_mc_global * sigma_y_mc_global + ...
                w_y_cs_global * sigma_y_cs_global;

            %% estimate disparity histogram for all image pairs by bootstrapping
            
            % calculate number of particles to be removed
            num_particles_remove = round(percentage_particles_remove * num_particles);
            fprintf('num particles remove: %d\n', num_particles_remove);

            % directory where results for this case will be saved
            current_results_directory = fullfile(top_results_directory, ['overlap=' logical_string{particle_overlap + 1}], ['delta_x=' num2str(delta_x, '%.1f') 'pix_std=' num2str(displacement_noise_std, '%.1f') 'pix_imnoise=' num2str(image_noise_std, '%02d')], [num2str(image_height) 'x' num2str(image_width)], ...
                ['ntot' num2str(num_particles, '%02d') '_npr' num2str(num_particles_remove, '%02d') '_nbs' num2str(num_bootstrap_trials)], ['particle_position_calc=' particle_position_calc_method]);
            if strcmp(particle_position_calc_method, 'id')
                current_results_directory = fullfile(current_results_directory, ['thresh=' num2str(intensity_threshold) '_size=' sizeprops.method]);
            end
            % create results directory
            mkdir_c(current_results_directory);

            fprintf('Estimating disparity by bootstrapping (local)\n');

            %% initialize variables
            w_x_im_local = nans(1, num_image_pairs);
            w_y_im_local = nans(1, num_image_pairs);

            w_x_mc_local = nans(1, num_image_pairs);
            w_y_mc_local = nans(1, num_image_pairs);

            w_x_cs_local = nans(1, num_image_pairs);
            w_y_cs_local = nans(1, num_image_pairs);

            sigma_x_im_local_avg = nans(1, num_image_pairs);
            sigma_y_im_local_avg = nans(1, num_image_pairs);
            sigma_x_mc_local_avg = nans(1, num_image_pairs);
            sigma_y_mc_local_avg = nans(1, num_image_pairs);
            sigma_x_cs_local_avg = nans(1, num_image_pairs);
            sigma_y_cs_local_avg = nans(1, num_image_pairs);

            sigma_x_weighted_local = nans(1, num_image_pairs);
            sigma_y_weighted_local = nans(1, num_image_pairs);

            %% loop through image pairs
            for image_pair_index = 1:num_image_pairs
                fprintf('image pair: %d\n', image_pair_index);
                % extract current image
                im1 = double(image_matrix_1(:, :, image_pair_index));
                im2 = double(image_matrix_2(:, :, image_pair_index));

                if strcmp(particle_position_calc_method, 'true')
                    % extract particle positions for the first image
                    x1 = particle_positions.ParticlePositions(image_pair_index).X1;
                    y1 = particle_positions.ParticlePositions(image_pair_index).Y1;

                    % extract particle positions for the second image
                    x2 = particle_positions.ParticlePositions(image_pair_index).X2;
                    y2 = particle_positions.ParticlePositions(image_pair_index).Y2;
                elseif strcmp(particle_position_calc_method, 'id')
                    %%
                    % calculate product of images
                    im_p = im1.*im2/max(im1(:));
                    
                    % identify intensity peaks using particle identification
                    [p_matrix,peaks,num_p] = dynamic_threshold_segmentation_v3(im_p, intensity_threshold, 0);
                    
                    % perform particle sizing
                    [XYDiameter, mapsizeinfo, locxy]=particle_size_MAIN_V1(im_p, p_matrix, num_p, sizeprops);
                    
                    % extract locations of particles
                    x1 = XYDiameter(:, 1);
                    y1 = XYDiameter(:, 2);
                    
                    % identify nan elements
                    indices = isnan(x1) | isnan(y1);
                    
                    % remove nan elements
                    x1(indices) = [];
                    y1(indices) = [];
                    
                    % assign the positions of the peaks in the second image
                    % to be the same (as the disparity is within subpixel
                    % limit)
                    x2 = x1;
                    y2 = y1;
                    %%
                else
                    fprintf('unknown option for particle position. exiting.\n');
                    return;
                end

                %% calculate local variance of the uncertainty schemes by bootstrapping
                % calculate variance of the uncertaintys chemes
                [sigma_avg, sigma_std] = calculate_variance_boostrapping(im1, im2, X, Y, U, V, x1, y1, x2, y2, d_p, num_particles_remove, num_bootstrap_trials, uncertainty_flags);
                
                % extract average of uncertainty schemes
                sigma_x_im_local_avg = sigma_avg.imx;
                sigma_y_im_local_avg = sigma_avg.imy;
                sigma_x_mc_local_avg = sigma_avg.mcx;
                sigma_y_mc_local_avg = sigma_avg.mcy;
                sigma_x_cs_local_avg = sigma_avg.csx;
                sigma_y_cs_local_avg = sigma_avg.csy;

                % extract variance of uncertainty schemes
                std_sigma_x_im_local = sigma_std.imx;
                std_sigma_y_im_local = sigma_std.imy;
                std_sigma_x_mc_local = sigma_std.mcx;
                std_sigma_y_mc_local = sigma_std.mcy;
                std_sigma_x_cs_local = sigma_std.csx;
                std_sigma_y_cs_local = sigma_std.csy;
                
                %% calculate local weights

                w_x_im_local(image_pair_index) = 1/std_sigma_x_im_local^2/(1/std_sigma_x_im_local^2 + 1/std_sigma_x_mc_local^2 + 1/std_sigma_x_cs_local^2);
                w_y_im_local(image_pair_index) = 1/std_sigma_y_im_local^2/(1/std_sigma_y_im_local^2 + 1/std_sigma_y_mc_local^2 + 1/std_sigma_y_cs_local^2);

                w_x_mc_local(image_pair_index) = 1/std_sigma_x_mc_local^2/(1/std_sigma_x_im_local^2 + 1/std_sigma_x_mc_local^2 + 1/std_sigma_x_cs_local^2);
                w_y_mc_local(image_pair_index) = 1/std_sigma_y_mc_local^2/(1/std_sigma_y_im_local^2 + 1/std_sigma_y_mc_local^2 + 1/std_sigma_y_cs_local^2);

                w_x_cs_local(image_pair_index) = 1 - (w_x_im_local(image_pair_index) + w_x_mc_local(image_pair_index));
                w_y_cs_local(image_pair_index) = 1 - (w_y_im_local(image_pair_index) + w_y_mc_local(image_pair_index));

                %% calculate weighted uncertainty
                
                % calculate weighted uncertainty for this image pair
                sigma_x_weighted_local(image_pair_index) = w_x_im_local(image_pair_index) * sigma_x_im_global(image_pair_index) + ...
                    w_x_mc_local(image_pair_index) * sigma_x_mc_global(image_pair_index) + ...
                    w_x_cs_local(image_pair_index) * sigma_x_cs_global(image_pair_index);

                sigma_y_weighted_local(image_pair_index) = w_y_im_local(image_pair_index) * sigma_y_im_global(image_pair_index) + ...
                    w_y_mc_local(image_pair_index) * sigma_y_mc_global(image_pair_index) + ...
                    w_y_cs_local(image_pair_index) * sigma_y_cs_global(image_pair_index);

            end

            %% remove invalid measurements

            % image matching
            sigma_im_global_all_valid = remove_invalid_measurements([sigma_x_im_global, sigma_y_im_global], min_uncertainty_threshold, max_uncertainty_threshold);

            % moment of correlation
            sigma_mc_global_all_valid = remove_invalid_measurements([sigma_x_mc_global, sigma_y_mc_global], min_uncertainty_threshold, max_uncertainty_threshold);

            % correlation statistics
            sigma_cs_global_all_valid = remove_invalid_measurements([sigma_x_cs_global, sigma_y_cs_global], min_uncertainty_threshold, max_uncertainty_threshold);

            % unweighted uncertainty
            sigma_weighted_none_all_valid = remove_invalid_measurements([sigma_x_weighted_none, sigma_y_weighted_none], min_uncertainty_threshold, max_uncertainty_threshold);

            % globally weighted uncertainty
            sigma_weighted_global_all_valid = remove_invalid_measurements([sigma_x_weighted_global, sigma_y_weighted_global], min_uncertainty_threshold, max_uncertainty_threshold);

            % locally weighted uncertainty
            sigma_weighted_local_all_valid = remove_invalid_measurements([sigma_x_weighted_local, sigma_y_weighted_local], min_uncertainty_threshold, max_uncertainty_threshold);

            %% calculate pdfs of uncertainties

            % calculate pdf of global uncertainty for image matching
            N_im_global = histcounts(sigma_im_global_all_valid, bins, 'normalization', 'pdf');

            % calculate pdf of global uncertainty for moment of correlation
            N_mc_global = histcounts(sigma_mc_global_all_valid, bins, 'normalization', 'pdf');

            % calculate pdf of global uncertainty for correlation statistics
            N_cs_global = histcounts(sigma_cs_global_all_valid, bins, 'normalization', 'pdf');

            % calculate pdf of unweighted uncertainty
            N_weighted_none = histcounts(sigma_weighted_none_all_valid, bins, 'normalization', 'pdf');

            % calculate pdf of globally weighted uncertainty
            N_weighted_global = histcounts(sigma_weighted_global_all_valid, bins, 'normalization', 'pdf');

            % calculate pdf of locally weighted uncertainty
            N_weighted_local = histcounts(sigma_weighted_local_all_valid, bins, 'normalization', 'pdf');

            %% calculate rms of uncertainties

            % calculate rms of true uncertainty
            sigma_rms_true = displacement_noise_std/sqrt(num_particles);

            % calculate rms of global uncertainty for image matching
            sigma_rms_im_global = nanrms(sigma_im_global_all_valid);

            % calculate rms of global uncertainty for moment of correlation
            sigma_rms_mc_global = nanrms(sigma_mc_global_all_valid);

            % calculate rms of global uncertainty for correlation statistics
            sigma_rms_cs_global = nanrms(sigma_cs_global_all_valid);

            % calculate rms of unweighted uncertainty
            sigma_rms_weighted_none = nanrms(sigma_weighted_none_all_valid);

            % calculate rms of globally weighted uncertainty
            sigma_rms_weighted_global = nanrms(sigma_weighted_global_all_valid);

            % calculate rms of locally weighted uncertainty
            sigma_rms_weighted_local = nanrms(sigma_weighted_local_all_valid);

            %% save results

            % save uncertainties from individual schemes
            save(fullfile(current_results_directory, 'sigma_individual.mat'), ...
                'sigma_im_global_all_valid', 'sigma_mc_global_all_valid', ...
                'sigma_cs_global_all_valid');

            % save weighted uncertainties
            save(fullfile(current_results_directory, 'sigma_combined.mat'), ...
                'sigma_weighted_none_all_valid', ...
                'sigma_weighted_global_all_valid', ...
                'sigma_weighted_local_all_valid');

            % save weights
            save(fullfile(current_results_directory, 'weights.mat'), ...
                'w_x_im_local', 'w_y_im_local', 'w_x_mc_local', 'w_y_mc_local', 'w_x_cs_local', 'w_y_cs_local', ...
                'w_x_im_global', 'w_y_im_global', 'w_x_mc_global', 'w_y_mc_global', 'w_x_cs_global', 'w_y_cs_local');

            % save pdfs
            save(fullfile(current_results_directory, 'histograms.mat'), ...
                'bins', 'N_im_global', 'sigma_rms_im_global', ...
                'N_mc_global', 'sigma_rms_mc_global', ...
                'N_cs_global', 'sigma_rms_cs_global', ...
                'N_weighted_none', 'sigma_rms_weighted_none', ...
                'N_weighted_global', 'sigma_rms_weighted_global', ...
                'N_weighted_local', 'sigma_rms_weighted_local', ...
                'sigma_rms_true');

            %% plot histogram of uncertainties

            % directory to save figures
            figure_save_directory = fullfile(current_results_directory, 'figures');
            mkdir_c(figure_save_directory);


            y_max = 100;

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

            legend([l_im, l_mc, l_cs, l_wt_none, l_wt_global, l_wt_local, l_true], ...
                {'\sigma_{IM}', '\sigma_{MC}', '\sigma_{CS}', 'Unweighted', 'Global', 'Local', 'RMS, True'}, 'fontsize', 12)

            xlabel('Uncertainty (pix.)')
            ylabel('Probability Density', 'fontsize', 12)

            set(gcf, 'Position', [360   573   467   327])
            save_figure_to_png_eps_fig(figure_save_directory, 'histograms-comparison', [1, 0, 0]);
            %%
            close all
        end
    end
end
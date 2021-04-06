% function [sigma_trials, snr_metric_trials] = calculate_uncertainty_resampling(im1, im2, x1, y1, x2, y2, d_p, num_particles_resample, num_resampling_trials, uncertainty_methods, uncertainty_flags, fill_intensity, window_size, window_resolution, display_resampled_images)
function [sigma_trials, snr_metric_trials, U_trials, V_trials] = calculate_uncertainty_resampling_general(im1, im2, d_p, sizeprops, resampling_method, ...
                            resampling_percentage, num_resampling_trials, uncertainty_methods, fill_intensity, current_pass_settings, display_resampled_images)
% Function to resample images and calculate uncertainty.
%
% INPUTS:
% im1, im2: deformed image pair
% d_p: diameter of a particle in the image
% num_particles_resample: number of particles to remove in each trial
% resampling_method: 'remove-paired' or 'remove-random' or 'add-random'
% num_resampling_trials: number of resampling trials
% uncertainty_methods: array of uncertainty method names
% fill_intensity: intensity value to fill the region containing the removed
% particle
% current_pass_settings: processing settings for the current pass
% display_resampled_images: display resampled images? (true/false)
%
% OUTPUTS:
% sigma_trials: uncertainties for each scheme for all trials
% snr_metric_trials: snr metrics for all trials
% U_trials, V_trials: displacements for all trials
%
% AUTHOR:
% Lalit Rajendran (lrajendr@purdue.edu)
    
    % number of uncertainty methods
    num_uncertainty_methods = numel(uncertainty_methods);

    % ================================================
    %% identify particles
    % ================================================            
    if strcmpi(resampling_method, 'remove-paired')
        % identify paired particles
        im = sqrt(im1 .* im2);
        [x, y, d, Imax, num_particles] = identify_particles_02(im, d_p, sizeprops);
        x1 = x; x2 = x;
        y1 = y; y2 = y;
        d1 = d; d2 = d;
        Imax1 = Imax; Imax2 = Imax;        
    else
        % identify all particles
        [x1, y1, d1, Imax1, num_particles_1] = identify_particles_02(im1, d_p, sizeprops);
        [x2, y2, d2, Imax2, num_particles_2] = identify_particles_02(im2, d_p, sizeprops);

        % calculate number of particles identified
        % num_particles = round(0.5 * (num_particles_1 + num_particles_2));
        num_particles = min(num_particles_1, num_particles_2);            
    end

    % calculate number particles for resampling
    num_particles_resample = round(resampling_percentage * num_particles);

    % calculate image size
    [image_height, image_width] = size(im1);

    % ==========================
    % initialize arrays
    % ==========================        
    % displacements
    U_trials = nans(1, num_resampling_trials);
    V_trials = nans(1, num_resampling_trials);
    % uncertainties
    sigma_temp = cell(1, num_resampling_trials);
    % snr metrics
    snr_metric_temp = cell(1, num_resampling_trials);

    % ==========================
    % loop through resampling trials
    % ==========================
    parfor trial_index = 1:num_resampling_trials            

        % create temporary variables for the images
        im1_temp = im1;
        im2_temp = im2;
        
        % ==========================
        % create list of particle co-ordinates for addition or removal
        % ==========================
        if strcmpi(resampling_method, 'remove-paired')
            % remove paired particles
            [x1_r, y1_r, d1_r, Imax1_r] = resample_elements(x1, y1, d1, Imax1, num_particles_resample);
            x2_r = x1_r; y2_r = y1_r; d2_r = d1_r; Imax2_r = Imax1_r;

        elseif strcmpi(resampling_method, 'remove-random')
            % remove random particles
            [x1_r, y1_r, d1_r, Imax1_r] = resample_elements(x1, y1, d1, Imax1, num_particles_resample);
            [x2_r, y2_r, d2_r, Imax2_r] = resample_elements(x2, y2, d2, Imax2, num_particles_resample);

        else
            % add random particles
            x1_r = rand(1, num_particles_resample) * image_width;
            y1_r = rand(1, num_particles_resample) * image_height;
            d1_r = mean(d1, 'omitnan') * ones(1, num_particles_resample);
            Imax1_r = mean(Imax1, 'omitnan') * ones(1, num_particles_resample);

            x2_r = rand(1, num_particles_resample) * image_width;
            y2_r = rand(1, num_particles_resample) * image_height;
            d2_r = mean(d2, 'omitnan') * ones(1, num_particles_resample);
            Imax2_r = mean(Imax2, 'omitnan') * ones(1, num_particles_resample);
            
            % fprintf('%d, %d, %d, %d \n', numel(x1_r), numel(y1_r), numel(x2_r), numel(y2_r));
        end            
        
        % ==========================
        % resample particle intensities
        % ==========================                        
        if contains(resampling_method, 'remove')
            % remove particles
            im1_temp = remove_particles_from_image(im1, x1_r, y1_r, d1_r, prctile(im1(:), 10)); %, fill_intensity
            im2_temp = remove_particles_from_image(im2, x2_r, y2_r, d2_r, prctile(im2(:), 10)); %, fill_intensity                
        else
            % add particles
            im1_temp = im1_temp + generateParticleImage(image_height, image_width, x1_r, y1_r, d1_r, Imax1_r);
            im2_temp = im2_temp + generateParticleImage(image_height, image_width, x2_r, y2_r, d2_r, Imax2_r);
        end

        % display image and mark particles that were removed
        if display_resampled_images
            display_resampled_particle_images(im1, im2, im1_temp, im2_temp, x1, y1, x2, y2, x1_r, y1_r, x2_r, y2_r);
        end
        
        % calculate planar uncertainties            
        [sigma_temp{trial_index}, snr_metric_temp{trial_index}, U_trials(trial_index), V_trials(trial_index)] = calculate_planar_uncertainties(im1_temp, im2_temp, ...
                                                                                                                            current_pass_settings, uncertainty_methods);        
    end    

    % ==========================
    % extract results
    % ==========================
    sigma_trials = struct;
    snr_metric_trials = struct;
    for trial_index = 1:num_resampling_trials
        % uncertainties
        for method_index = 1:num_uncertainty_methods
            sigma_trials.([lower(uncertainty_methods{method_index}) 'x'])(trial_index) = sigma_temp{trial_index}{method_index}.x;
            sigma_trials.([lower(uncertainty_methods{method_index}) 'y'])(trial_index) = sigma_temp{trial_index}{method_index}.y;
        end        

        % SNR
        snr_metric_trials.PPR(trial_index) = snr_metric_temp{trial_index}.PPR;
        snr_metric_trials.MI(trial_index) = snr_metric_temp{trial_index}.MI;
    end
end

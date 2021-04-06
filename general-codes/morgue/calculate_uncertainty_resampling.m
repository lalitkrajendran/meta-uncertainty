% function [sigma_trials, snr_metric_trials] = calculate_uncertainty_resampling(im1, im2, x1, y1, x2, y2, d_p, num_particles_remove, num_resampling_trials, uncertainty_methods, uncertainty_flags, fill_intensity, window_size, window_resolution, display_resampled_images)
function [sigma_trials, snr_metric_trials, U_trials, V_trials] = calculate_uncertainty_resampling(im1, im2, x1, y1, x2, y2, d_p, num_particles_remove, num_resampling_trials, uncertainty_methods, fill_intensity, current_pass_settings, display_resampled_images)
% Function to calculate variance of an uncertainty schemes for a pair of
% deformed images with a known velocity field. The variance is calcualted 
% by randomly removing a few particles everytime and calculating the
% uncertainty by jack-knife resampling.
%
% INPUTS:
% im1, im2: deformed image pair
% x1, y1: location of particles in im1 (arrays, double)
% x2, y2: location of particles in im2 (arrays, double)
% d_p: diameter of a particle in the image
% num_particles_remove: number of particles to remove in each trial
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

    % initialize arrays
    sigma_x_local_trial = nans(num_uncertainty_methods, num_resampling_trials);
    sigma_y_local_trial = nans(num_uncertainty_methods, num_resampling_trials);
    ppr_trials = nans(1, num_resampling_trials);
    mi_trials = nans(1, num_resampling_trials);

    U_trials = nans(1, num_resampling_trials);
    V_trials = nans(1, num_resampling_trials);
    sigma_x_im_local_trial = nans(1, num_resampling_trials);
    sigma_y_im_local_trial = nans(1, num_resampling_trials);
    sigma_x_mc_local_trial = nans(1, num_resampling_trials);
    sigma_y_mc_local_trial = nans(1, num_resampling_trials);
    sigma_x_cs_local_trial = nans(1, num_resampling_trials);
    sigma_y_cs_local_trial = nans(1, num_resampling_trials);
    
    % calculate image size
    [image_height, image_width] = size(im1);

    % fill_intensity = prctile(im1(:), 25);
    %% loop through trials
    parfor trial_index = 1:num_resampling_trials
        % fprintf('resampling trial: %d\n');
        % randomly sample particles for removal
        particle_sample_1 = randsample(numel(x1), num_particles_remove);
        particle_sample_2 = randsample(numel(x2), num_particles_remove);
        
        %% remove the pixels surrounding these particles from the images
        im1_temp = im1;
        im2_temp = im2;

        for sample_index = 1:num_particles_remove

            % index of the particle that is to be removed
            particle_index_1 = particle_sample_1(sample_index);
            particle_index_2 = particle_sample_2(sample_index);

            % x, y location of the particle to be removed
            x1_current = x1(particle_index_1);
            y1_current = y1(particle_index_1);
            x2_current = x2(particle_index_2);
            y2_current = y2(particle_index_2);
            
            % diameter of the particle to be removed 
            d_current = round(d_p(particle_index_1));
            
            % integer locations on the image
            % r1_array = (-d_current:d_current) + image_height - round(y1_current) + 1;
            r1_array = (-d_current:d_current) + round(y1_current);
            c1_array = (-d_current:d_current) + round(x1_current);
            % r2_array = (-d_current:d_current) + image_height - round(y2_current) + 1;
            r2_array = (-d_current:d_current) + round(y2_current);
            c2_array = (-d_current:d_current) + round(x2_current);

            % loop through locations and set pixel intensities to zero
            for r_index = 1:numel(r1_array)
                for c_index = 1:numel(c1_array)

                    % row, col locations for the current particle in im1
                    r1_current = round(r1_array(r_index));
                    c1_current = round(c1_array(c_index));

                    % zero out intensity map
                    if check_valid_index(image_height, image_width, r1_current, c1_current)
                        im1_temp(r1_current, c1_current) = fill_intensity;
                    end

                    % row, col locations for the current particle in im2
                    r2_current = round(r2_array(r_index));
                    c2_current = round(c2_array(c_index));

                    % zero out intensity map
                    if check_valid_index(image_height, image_width, r2_current, c2_current)
                        im2_temp(r2_current, c2_current) = fill_intensity;
                    end
                end
            end
        end

        %% display image and mark particles that were removed
        if display_resampled_images
            clim = [0, prctile([im1(:); im2(:)], 90)]; %0.5 * max(im1(:))];
            
            figure
            subplot(2, 2, 1)
            imagesc(im1)
            colormap(gray)
            hold on
            plot(x1, y1, 'go')
            plot(x1(particle_sample_1), y1(particle_sample_1), 'ro')
            caxis(clim)
            colorbar
            set_axes(gca);
            axis off
            title('Frame 1, Original')
            
            subplot(2, 2, 2)
            imagesc(im1_temp)
            colormap(gray)
            hold on
            plot(x1, y1, 'go')
            plot(x1(particle_sample_1), y1(particle_sample_1), 'ro')
            caxis(clim)
            colorbar
            set_axes(gca);
            axis off
            title('Frame 1, Resampled')

            subplot(2, 2, 3)
            imagesc(im2)
            colormap(gray)
            hold on
            plot(x2, y2, 'go')
            plot(x2(particle_sample_2), y2(particle_sample_2), 'ro')
            caxis(clim)
            colorbar
            set_axes(gca);
            axis off
            title('Frame 2, Original')
            
            subplot(2, 2, 4)
            imagesc(im2_temp)
            colormap(gray)
            hold on
            plot(x2, y2, 'go')
            plot(x2(particle_sample_2), y2(particle_sample_2), 'ro')
            caxis(clim)
            colorbar
            set_axes(gca);
            axis off
            title('Frame 2, Original')
            
            set(gcf, 'Position', [57           1        1000         700]);        
        end

        % calculate planar uncertainties
        % [sigma_current, snr_metric_current] = calculate_planar_uncertainties(im1_temp, im2_temp, window_size, window_resolution, uncertainty_methods, uncertainty_flags);
        [sigma_current, snr_metric_current, U_trials(trial_index), V_trials(trial_index)] = calculate_planar_uncertainties(im1_temp, im2_temp, current_pass_settings, uncertainty_methods);
        
        % extract uncertainties
        for method_index = 1:num_uncertainty_methods
            sigma_x_local_trial(method_index, trial_index) = sigma_current{method_index}.x;
            sigma_y_local_trial(method_index, trial_index) = sigma_current{method_index}.y;
        end

        % extract SNR
        ppr_trials(trial_index) = snr_metric_current.PPR;
        mi_trials(trial_index) = snr_metric_current.MI;
    end
    
    % ================================================
    %% store uncertainties from trials in a structure
    % ================================================
    sigma_trials = struct;
    
    for method_index = 1:num_uncertainty_methods
        sigma_trials.([lower(uncertainty_methods{method_index}) 'x']) = sigma_x_local_trial(method_index, :);
        sigma_trials.([lower(uncertainty_methods{method_index}) 'y']) = sigma_y_local_trial(method_index, :);
    end
   
    % store snr metric
    snr_metric_trials = struct;
    snr_metric_trials.PPR = ppr_trials;
    snr_metric_trials.MI = mi_trials;

end

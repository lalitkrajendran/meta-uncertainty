function [sigma_trials, snr_metric_trials] = calculate_uncertainty_resampling(im1, im2, x1, y1, x2, y2, d_p, num_particles_remove, num_resampling_trials, uncertainty_flags, fill_intensity, window_size, window_resolution, display_resampled_images)
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
% uncertainty_flags: structure containing flags on which uncertainty
% schemes to evaluate. The MC flag needs to be turned on
% fill_intensity: intensity value to fill the region containing the removed
% particle
% window_size: window size for processing [x, y] (pix.)
% window_resolution: window resolution for processing [x, y] (pix.)
% display_resampled_images: display resampled images? (true/false)
%
% OUTPUTS:
% sigma_trials: uncertainties for each scheme for all trials
% snr_metric_trials: snr metrics for all trials
%
% AUTHOR:
% Lalit Rajendran (lrajendr@purdue.edu)

    % initialize arrays
    sigma_x_im_local_trial = nans(1, num_resampling_trials);
    sigma_y_im_local_trial = nans(1, num_resampling_trials);
    sigma_x_mc_local_trial = nans(1, num_resampling_trials);
    sigma_y_mc_local_trial = nans(1, num_resampling_trials);
    sigma_x_cs_local_trial = nans(1, num_resampling_trials);
    sigma_y_cs_local_trial = nans(1, num_resampling_trials);
    ppr_trials = nans(1, num_resampling_trials);
    mi_trials = nans(1, num_resampling_trials);
    
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

        % ================================================
        %% calculate planar uncertainties
        % ================================================
        % calculate uncertainty using IM        
        % [sigma_x_im_local_trial(trial_index), sigma_y_im_local_trial(trial_index), ~, ~] = original_particle_disparity_no_deform(im1_temp, im2_temp, image_width/2, image_height/2, 0, 0, image_height);
        [sigma_x_im_local_trial(trial_index), sigma_y_im_local_trial(trial_index), ~, ~] = original_particle_disparity(im1_temp, im2_temp, image_width/2, image_height/2, 0, 0, window_resolution(1), 0);

        % calculate uncertainty using MC        
        [~, ~, U, V, ~, ~, ~, uncertainty2D, SNRmetric] = PIVwindowed(im1_temp, im2_temp, 'SCC', window_size, [window_resolution; window_resolution], 0, [2.8, 2.8], 1, 3, 0, 0, 0, image_width/2, image_height/2, uncertainty_flags, 1, 0, 0);
        % [~,~,Uc,Vc,~,~,~,uncertainty2D,~]=PIVwindowed(im1d,im2d,Corr{e},Wsize(e,:),Wres(:, :, e),0,D(e,:),Zeromean(e),Peaklocator(e),find_extrapeaks,frac_filt(e),saveplane(e),X(Eval>=0),Y(Eval>=0),uncertainty,e);
        
        sigma_x_mc_local_trial(trial_index) = sqrt(uncertainty2D.biasx.^2 + (uncertainty2D.Ixx.^2) ./ uncertainty2D.Neff);
        sigma_y_mc_local_trial(trial_index) = sqrt(uncertainty2D.biasy.^2 + (uncertainty2D.Iyy.^2) ./ uncertainty2D.Neff);

        % extract snr metric
        ppr_trials(trial_index) = SNRmetric.PPR;
        mi_trials(trial_index) = SNRmetric.MI;

        % calculate uncertainty using CS
        [sigma_x_cs_local_trial(trial_index), sigma_y_cs_local_trial(trial_index)] = correlation_statistics(im1_temp, im2_temp, window_size, [window_resolution; window_resolution], window_size(1)/2, window_size(2)/2);        
    end
    
    % ================================================
    %% store uncertainties from trials in a structure
    % ================================================
    sigma_trials = struct;
    
    sigma_trials.imx = sigma_x_im_local_trial;
    sigma_trials.imy = sigma_y_im_local_trial;
    
    sigma_trials.mcx = sigma_x_mc_local_trial;
    sigma_trials.mcy = sigma_y_mc_local_trial;
    
    sigma_trials.csx = sigma_x_cs_local_trial;
    sigma_trials.csy = sigma_y_cs_local_trial;
   
    % store snr metric
    snr_metric_trials = struct;
    snr_metric_trials.ppr = ppr_trials;
    snr_metric_trials.mi = mi_trials;
end

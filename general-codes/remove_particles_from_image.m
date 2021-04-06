function im_new = remove_particles_from_image(im, xp, yp, dp, fill_intensity)
% Function to remove particle intensity maps from an image and replace that
% region with a pre-defined intensity level
%
% INPUTS:
% im: image where the particles are to be removed
% xp, yp: centroids of the particles to be remoed (pix.)
% dp: diameters of the particles to be removed (pix.)
% fill_intensity: intensity to be assigned to regions where particles will be removed
%
% OUTPUTS:
% im_new: image after particle removal
%
% AUTHOR:
% Lalit Rajendran (lrajendr@purdue.edu)

    % get image size
    [image_height, image_width] = size(im);

    % calculate number of particles
    num_particles = numel(xp);

    % copy image intensity to new array
    im_new = im;

    % loop through particles and remove intensity maps
    for particle_index = 1:num_particles        
        % generate array of particle extent        
        [r_array, c_array] = generate_particle_extent_array(xp(particle_index), yp(particle_index), dp(particle_index));

        % loop through locations and set pixel intensities to fill intensity
        for r_index = 1:numel(r_array)
            for c_index = 1:numel(c_array)
                % row, col locations for the current particle in im1
                r_current = round(r_array(r_index));
                c_current = round(c_array(c_index));

                % replace intensity map with fill intensity
                if check_valid_index(image_height, image_width, r_current, c_current)
                    im_new(r_current, c_current) = fill_intensity;
                end
            end
        end
    end
end
function [unc_planar, snr_metric] = calculate_planar_uncertainties(im1, im2, window_size, window_resolution, uncertainty_methods, uncertainty_flags)

    % number of uncertainty methods
    num_methods = numel(uncertainty_methods);

    % calculate image size
    [image_height, image_width] = size(im1);

    % create cell to hold uncertainties
    unc_planar = cell(1, num_methods);

    for method_index = 1:num_methods
        unc_planar{method_index} = struct;
        unc_planar{method_index}.name = uncertainty_methods{method_index};
            
        % calculate uncertainty using IM        
        if strcmp(uncertainty_methods{method_index}, 'IM')
            [sigma_x, sigma_y, ~, ~] = original_particle_disparity(im1, im2, image_width/2, image_height/2, 0, 0, window_resolution(1), 0);
        % calculate uncertainty using MC
        elseif strcmp(uncertainty_methods{method_index}, 'MC')
            % [~, ~, U, V, ~, ~, ~, uncertainty2D, SNRmetric] = PIVwindowed(im1, im2, 'SCC', window_size, [window_resolution; window_resolution], 0, [2.8, 2.8], 1, 3, 0, 0, 0, image_width/2, image_height/2, uncertainty_flags, 1, 0, 0);            
            [~, ~, U, V, ~, ~, ~, uncertainty2D, SNRmetric] = PIVwindowed(im1, im2, 'SCC', window_size, [window_resolution; window_resolution], 0, [2.8, 2.8], 1, 1, 0, 1, 0, image_width/2, image_height/2, uncertainty_flags, 1, 0, 0);            
            sigma_x = sqrt(uncertainty2D.biasx.^2 + (uncertainty2D.Ixx.^2) ./ uncertainty2D.Neff);
            sigma_y = sqrt(uncertainty2D.biasy.^2 + (uncertainty2D.Iyy.^2) ./ uncertainty2D.Neff);
            snr_metric = SNRmetric;
        % calculate uncertainty using CS
        elseif strcmp(uncertainty_methods{method_index}, 'CS')
            [sigma_x, sigma_y] = correlation_statistics(im1, im2, window_size, [window_resolution; window_resolution], window_size(1)/2, window_size(2)/2);        
        end

        % allocate uncertainties
        unc_planar{method_index}.x = sigma_x;
        unc_planar{method_index}.y = sigma_y;
    end
    
end

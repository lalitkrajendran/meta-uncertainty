function [x1, y1, x2, y2] = identify_particles(im1, im2, d_p, sizeprops, display_id_results)
    % ================================================
    %% identify particles
    % ================================================                
    % calculate product of images
    im_p = sqrt(im1 .* im2); %/sqrt(max(im1(:)) * max(im2(:)));
    
    % set minimum intensity threshold for the image    
    intensity_threshold = prctile(im_p(:), 90);

    % identify intensity peaks using particle identification
    [p_matrix, peaks, num_p] = dynamic_threshold_segmentation_v3(im_p, intensity_threshold, 0);

    % ================================================
    %% dot sizing and centroid estimation
    % ================================================                
    % change size properties
    % particle sizing settings
    sizeprops_current = sizeprops;
    sizeprops_current.p_area = d_p; %0.5*d_p;
    
    % perform particle sizing
    [XYDiameter, mapsizeinfo, locxy] = particle_size_MAIN_V1(im_p, p_matrix, num_p, sizeprops_current);
    
    % extract locations of particles
    x1 = XYDiameter(:, 1);
    y1 = XYDiameter(:, 2);
    
    % extract particle diameters
    d1 = XYDiameter(:, 3);
    
    % identify nan elements
    indices = isnan(x1) | isnan(y1);
    
    % remove nan elements
    x1(indices) = [];
    y1(indices) = [];
    d1(indices) = [];
    
    % assign the positions of the peaks in the second image
    % to be the same (as the disparity is within subpixel
    % limit)
    x2 = x1;
    y2 = y1;
    
    if display_id_results
        % ================================================
        %% display identification results    
        % ================================================
        figure
        imagesc(im_p)
        colormap(gray)
        caxis([0 intensity_threshold*2])
        % maxval = 0.8 * max(im_p(:))/2;
        % caxis([0 maxval])
        hold on
        plot(x1, y1, 'o')
        
        % plot window edge
        x_c = size(im_p, 2)/2;
        y_c = size(im_p, 1)/2;
        
        x = [x_c - x_c/2, x_c + x_c/2, x_c + x_c/2, x_c - x_c/2, x_c - x_c/2];
        y = [y_c - y_c/2, y_c - y_c/2, y_c + y_c/2, y_c + y_c/2, y_c - y_c/2];
        plot(x, y, '*-')
        colorbar
        set_axes(gca);
        pause(0.1);            
    end

end
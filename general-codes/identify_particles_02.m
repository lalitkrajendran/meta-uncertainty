function [x, y, d, Imax, N] = identify_particles_02(im_p, d_p, sizeprops, display_id_results)
% Function to identify and size particles on an interrogation window
%
% INPUTS: 
% im_p: image
% d_p: expected particle diameter
% sizeprops: sizing settings
% display_id_results: True or False
%
% OUTPUTS:
% x, y: particle locations (pix.)
% d: particle diameter (pix.)
% N: number of particles
% 
% AUTHORS:
% Lalit Rajendran (lrajendr@purdue.edu)
%
% DATE:
% 2020/06/05

    if nargin < 4
        display_id_results = false;
    end
    
    % ================================================
    %% identify particles
    % ================================================                
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
    x = XYDiameter(:, 1);
    y = XYDiameter(:, 2);
    
    % extract particle diameters
    d = XYDiameter(:, 3);
    
    % extract max intensities
    Imax = XYDiameter(:, 4);

    % % identify nan elements
    % indices = isnan(x) | isnan(y);
    
    % % remove nan elements
    % x(indices) = [];
    % y(indices) = [];
    % d(indices) = [];
    % Imax(indices) = [];

    [x, y, d, Imax] = remove_nan_elements(x, y, d, Imax);

    % number of particles
    N = numel(x);
    
    % ================================================
    %% display identification results    
    % ================================================
    if display_id_results        
        figure
        imagesc(im_p)
        colormap(gray)
        caxis([0 intensity_threshold*2])
        % maxval = 0.8 * max(im_p(:))/2;
        % caxis([0 maxval])
        hold on
        plot(x, y, 'o')
        
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
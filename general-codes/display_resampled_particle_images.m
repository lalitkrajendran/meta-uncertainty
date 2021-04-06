function display_resampled_particle_images(im1, im2, im1_temp, im2_temp, x1, y1, x2, y2, x1_r, y1_r, x2_r, y2_r)
% Function to display images of particles that were resampled in the two frames
%
% INPUTS:
% im1, im2: original images 
% im1_temp, im2_temp: resampled images
% x1, y1, x2, y2: particle centroids on the original images 
% x1_r, y1_r, x2_r, y2_r: resampled particle centroids
% 
% OUTPUTS:
% None
%
% AUTHOR:
% Lalit Rajendran (lrajendr@purdue.edu)

    
    colors = lines(7);

    % c1 = colors(7, :);
    % c1 = colors(2, :);
    c1 = [1, 0, 0];
    c2 = colors(5, :);
    c3 = colors(1, :);

    % % calculate limits for colobar
    % clim = [0, prctile([im1(:); im2(:)], 90)]; %0.5 * max(im1(:))];
    % % axis limits
    % [image_height, image_width] = size(im1);
    % x_lim = [0.25 * image_width, 0.75 * image_width];
    % y_lim = [0.25 * image_height, 0.75 * image_height];

    % figure    
    % % Frame 1, Original
    % subplot(2, 2, 1)
    % imagesc(im1)
    % colormap(flipud(gray))
    % hold on
    % plot(x1, y1, 'o', 'color', c2)
    % plot(x1_r, y1_r, 'o', 'color', c1)
    % caxis(clim)
    % % colorbar
    % set_axes(gca);
    % axis off
    % title('Frame 1, Original')
    
    % % Frame 1, Resampled
    % subplot(2, 2, 2)
    % imagesc(im1_temp)
    % colormap(flipud(gray))
    % hold on
    % plot(x1, y1, 'o', 'color', c2)
    % plot(x1_r, y1_r, 'o', 'color', c1)
    % caxis(clim)
    % % colorbar
    % set_axes(gca);
    % axis off
    % title('Frame 1, Resampled')

    % % Frame 2, Original
    % subplot(2, 2, 3)
    % imagesc(im2)
    % colormap(flipud(gray))
    % hold on
    % plot(x2, y2, 'o', 'color', c2)
    % plot(x2_r, y2_r, 'o', 'color', c1)
    % caxis(clim)
    % % colorbar
    % set_axes(gca);
    % axis off
    % title('Frame 2, Original')
    
    % % Frame 2, Resampled
    % subplot(2, 2, 4)
    % imagesc(im2_temp)
    % colormap(flipud(gray))
    % hold on
    % plot(x2, y2, 'o', 'color', c2)
    % plot(x2_r, y2_r, 'o', 'color', c1)
    % caxis(clim)
    % % colorbar
    % set_axes(gca);
    % axis off
    % title('Frame 2, Original')
    
    % % set figure position
    % set(gcf, 'resize', 'off')
    % % set(gcf, 'Position', [57           1        1000         700]);           
    % set(gcf, 'Position', [114   290   557   492]);           
    % set(gcf, 'resize', 'off')
    % drawnow();

    % axis limits
    clim = [0, 100];
    [image_height, image_width] = size(im1);
    % x_lim = [0.25 * image_width, 0.75 * image_width];
    % y_lim = [0.25 * image_height, 0.75 * image_height];    
    x_lim = [0, image_width];
    y_lim = [0, image_height];    
    
    figure    
    for frame_index = 1:2
        if frame_index == 1
            im = im1;
            im_temp = im1_temp;
            x = x1_r;
            y = y1_r;
        else
            im = im2;
            im_temp = im2_temp;
            x = x2_r;
            y = y2_r;
        end

        % Frame 1, Original
        subplot(2, 2, (frame_index - 1) * 2 + 1)
        imagesc(im)
        colormap(flipud(gray))
        hold on
        % plot(x1, y1, 'o', 'color', c2)
        plot(x, y, 'o', 'color', c1, 'linewidth', 1.5)
        caxis(clim)
        % colorbar
        set_axes(gca);
        axis([x_lim, y_lim])
        axis off
        title(['Frame ' num2str(frame_index) ', Original'])
        
        % Frame 1, Resampled
        subplot(2, 2, (frame_index - 1) * 2 + 2)
        imagesc(im_temp)
        colormap(flipud(gray))        
        hold on
        % plot(x1, y1, 'o', 'color', c2)
        plot(x, y, 'o', 'color', c1, 'linewidth', 1.5)
        caxis(clim)
        % colorbar
        set_axes(gca);
        axis([x_lim, y_lim])
        axis off
        title(['Frame ' num2str(frame_index) ', Perturbed'])
    end

    % set figure position
    set(gcf, 'resize', 'off')
    % set(gcf, 'Position', [57           1        1000         700]);           
    % set(gcf, 'units', 'inches', 'Position', [327   455   535   270]/113);           
    set(gcf, 'units', 'inches', 'Position', [327   455   535   270*2]/113);           
    set(gcf, 'resize', 'off')
    drawnow();
end
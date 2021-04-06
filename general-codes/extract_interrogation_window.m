function zone = extract_interrogation_window(im, X, Y, window_size)
% Function to extract an interrogation window from an image.
%
% INPUTS:
% im - image
% X, Y - co-ordinates of the center of the image window
% window_size - size of the window to be extracted
%
% OUTPUTS:
% zone: interrogation window that is zero padded appropriately
%
% AUTHORS:
% Lalit Rajendran (lrajendr@purdue.edu)
%
% DATE:
% 06/12/2020

    % calculate image size
    [image_height, image_width] = size(im);
    
    % extract window size
    window_size_x = window_size(1);
    window_size_y = window_size(2);
    
    % identify window extents
    xmin = X - window_size_x/2 + 1;
    xmax = X + window_size_x/2;
    ymin = Y - window_size_y/2 + 1;
    ymax = Y + window_size_y/2;
    
    % find the image windows
    zone = im(max([1 ymin]):min([image_height ymax]), max([1 xmin]):min([image_width xmax]));
    
    % ================================================
    % zero pad windows
    % ================================================
    if size(zone,1)~=window_size_y || size(zone,2)~=window_size_x
        w = zeros(window_size_y,window_size_x);
        w( 1+max([0 1-ymin]):window_size_y - max([0 ymax-image_height]), 1+max([0 1-xmin]):window_size_x - max([0 xmax-image_width]) ) = zone;
        zone = w;
    end
    
end
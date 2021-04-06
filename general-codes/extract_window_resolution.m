function [window_resolution_x, window_resolution_y] = extract_window_resolution(pass_settings)
% Function to extract window resolutoin from a jobfile structure
%
% INPUT:
% pass_settings: jobfile parameters for this pass
%
% OUTPUTS:
% window_resolution_x, window_resolution_y: window resolutions along the x and y directions (pix.)
%
% AUTHOR:
% Lalit Rajendran (lrajendr@purdue.edu)
%
% DATE:
% 2020/06/12

    % split string
    str = strsplit(pass_settings.winres, ';');
    str2 = strsplit(str{1}, ',');
    % extract window resolution along x
    window_resolution_x = str2double(str2{1});
    % extract window resolution along y
    window_resolution_y = str2double(str2{2});
end

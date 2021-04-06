function [window_size_x, window_size_y] = extract_window_size(pass_settings)
% Function to extract window size from a jobfile structure
%
% INPUT:
% pass_settings: jobfile parameters for this pass
%
% OUTPUTS:
% window_size_x, window_size_y: window sizes along the x and y directions (pix.)
%
% AUTHOR:
% Lalit Rajendran (lrajendr@purdue.edu)
%
% DATE:
% 2020/06/12

    % split the string at the comma
    str = strsplit(pass_settings.winsize, ',');
    % extract window size along x
    window_size_x = str2double(str{1});
    % extract window size along y
    window_size_y = str2double(str{2});
end

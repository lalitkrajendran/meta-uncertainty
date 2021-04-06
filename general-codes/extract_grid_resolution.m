function [grid_resolution_x, grid_resolution_y] = extract_grid_resolution(pass_settings)
% Function to extract grid resolution from a jobfile structure
%
% INPUT:
% pass_settings: jobfile parameters for this pass
%
% OUTPUTS:
% grid_resolution_x, grid_resolution_y:grid resolutions along the x and y directions (pix.)
%
% AUTHOR:
% Lalit Rajendran (lrajendr@purdue.edu)
%
% DATE:
% 2020/06/18

    % split string
    str = strsplit(pass_settings.gridres, ',');
    % extract window resolution along x
    grid_resolution_x = str2double(str{1});
    % extract window resolution along y
    grid_resolution_y = str2double(str{2});

end
function [d_min, ind1, ind2] = calculate_minimum_distance_between_two_curves(x1, y1, x2, y2)
% Function to calculate minimum distance between two curves
%
% INPUTS:
% x1, y1: points on the first curve
% x2, y2: points on the second curve
%
% OUTPUTS:
% d_min: minimum distance
% ind1, ind2: indices of the pair of poins with the minimum distance
%
% AUTHORS:
% Lalit Rajendran (lrajendr@purdue.edu)
%
% DATE:
% 05/01/2020

    % calculate distance between every pair of points
    d = sqrt((x1 - x2.').^2 + (y1 - y2.').^2);
    % find minimum distance entry
    [d_min, indmin] = min(d(:));
    % calculate corresponding indices on the two arrays
    [ind1, ind2] = ind2sub([size(d, 1), size(d, 2)], indmin);
end
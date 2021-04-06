function [e, p] = calculate_joint_entropy(x1, x2, bins1, bins2)
% This function calculates the joint shannon entropy of a vector of observations.
% 
% INPUTS:
% x1, x2: vector of observations
% bins1, bins2: bins for the histogram calculation. if it is empty, then the
% default binning option is used
%
% OUTPUTS:
% e: entropy
% p: joint pdf
%
% AUTHOR:
% Lalit Rajendran (lrajendr@purdue.edu)

    if nargin < 3
        bins1 = [];
        bins2 = [];
    end
    
    % only retain real measurements
    x1 = real(x1);
    x2 = real(x2);

    % remove nan measurements
    r = find(isnan(x1) | isnan(x2));
    x1(r) = [];
    x2(r) = [];

    % calculate joint pdf
    if isempty(bins1) || isempty(bins2)
        [p, bins1, bins2] = histcounts2(x1, x2, 'normalization', 'pdf');
    else
        p = histcounts2(x1, x2, 'xbinedges', bins1, 'ybinedges', bins2, 'normalization', 'pdf');
    end
    
    % calculate shannon entropy
    h = p .* log(p);
    h(isnan(h)) = 0;
    e = -trapz(bins1(1:end-1), trapz(bins2(1:end-1), h));
end
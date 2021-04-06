function [e, p] = calculate_self_entropy(x, bins)
% This function calculates the shannon entropy of a vector of observations.
% 
% INPUTS:
% x: vector of observations
% bins: bins for the histogram calculation. if it is empty, then the
% default binning option is used
%
% OUTPUTS:
% e: entropy
% p: pdf
%
% AUTHOR:
% Lalit Rajendran (lrajendr@purdue.edu)

    if nargin < 2
        bins = [];
    end
    
    % only retain real measurements
    x = real(x);

    % remove nan measurements
    x(isnan(x)) = [];
    
    % calculate pdfs
    if isempty(bins)
        [p, bins] = histcounts(x, 'normalization', 'pdf');
    else
        p = histcounts(x, bins, 'normalization', 'pdf');
    end
    
    % calculate shannon entropy
    h = p .* log(p);
    h(isnan(h)) = 0;
    e = -trapz(bins(1:end-1), h);
end
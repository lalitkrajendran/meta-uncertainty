function [e, p] = calculate_shannon_entropy(x, bins)
    % This function calculates the shannon entropy of a vector of observations,
    % given by the formula: e = -sum(p * log2(p)), where p is the vector of
    % probabilities and e is the entropy
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
            % p = histcounts(x, 'normalization', 'probability');
            [p, bins] = histcounts(x, 'normalization', 'pdf');
        else
            % p = histcounts(x, bins, 'normalization', 'probability');
            p = histcounts(x, bins, 'normalization', 'pdf');
        end
        
        % calculate shannon entropy
        h = p .* log(p);
        h(isnan(h)) = 0;
        % e = -sum(p .* log2(p), 'omitnan');
        e = -trapz(bins(1:end-1), h);
    end
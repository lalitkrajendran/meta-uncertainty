function p_x0 = calculate_probability(x, x0, bin_edges)
% This function calculates the probability of a value x0 to lie in the
% set of values of x0
%
% INPUTS:
% x: array of observations
% x0: point at which probability is to be evaluated
% bin_edges: set of intervals
%
% OUTPUTS:
% p_x0: probability corresponding to this point

    % initialize bin edges
    if nargin < 3
        bin_edges = [];
    end
    
    % ensure values are real
    x = real(x);
    x0 = real(x0);
    
    % calculate probability distribution
    if isempty(bin_edges)
        [p, bin_edges] = histcounts(x, 'normalization', 'probability');
    else        
        [p, bin_edges] = histcounts(x, 'binedges', bin_edges, 'normalization', 'probability');
    end
    % calculate corresponding bin
    bin_index = find(x0 > bin_edges, 1, 'last');
    
    % calculate corresponding probability
    if ~isempty(bin_index) && bin_index < numel(p)
        p_x0 = p(bin_index);
    else
        p_x0 = 0;
    end
end
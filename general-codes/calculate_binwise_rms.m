function [indices, num_indices, err_rms_binned, unc_indiv_rms_binned, unc_comb_rms_binned] = calculate_binwise_rms(x, err, unc_indiv, unc_comb, bins);
% Function to calculate binwise rms errors and uncertainties corresponding to specified 
% ranges of the input variable
%
% INPUTS:
% x: variable that is to be used for binning
% err: error
% unc_indiv: individual uncertainties
% unc_comb: combined uncertainties
% bins: bin edges to use for discretization
% 
% OUTPUTS:
% indices: indices of elements in each bin
% num_indices: number of elements in each bin
% err_rms_binned: rms of errors in each bin
% unc_indiv_rms_binned, unc_comb_rms_binned: rms of individual and combined uncertainties in each bin
%
% AUTHOR:
% Lalit Rajendran (lrajendr@purdue.edu)
    
    % calculate number of individual methods
    num_individual_methods = size(unc_indiv, 2);
    % calculate number of combined methods
    num_combined_methods = size(unc_comb, 2);
    % calculate number of bins
    num_bins = numel(bins);
    
    % --------------------------
    % initialize variables
    % --------------------------
    indices = cell(num_bins, 1);
    num_indices = nans(num_bins, 1);
    err_rms_binned = nans(num_bins, 1);
    err_random_binned = nans(num_bins, 1);
    err_total_binned = nans(num_bins, 1);
    unc_indiv_rms_binned = nans(num_bins, num_individual_methods);
    unc_comb_rms_binned = nans(num_bins, num_individual_methods);

    % --------------------------
    % calculate binwise rms
    % --------------------------
    % discretize measurements into bins
    Y = discretize(x, bins);
    
    % loop through bins and calculate rms
    for bin_index = 1:num_bins
        % extract indices of elements that are in the current bin
        indices{bin_index} = find(Y == bin_index);
        % calculate number of elements in this bin
        num_indices(bin_index) = numel(indices{bin_index});
        % calculate rms error
        err_rms_binned(bin_index) = nanrms(err(indices{bin_index})');
        % calculate rms uncertainty - individual
        unc_indiv_rms_binned(bin_index, :) = nanrms(unc_indiv(indices{bin_index}, :), 1);
        % calculate rms uncertainty - combined
        unc_comb_rms_binned(bin_index, :) = nanrms(unc_comb(indices{bin_index}, :), 1);
    end    
end

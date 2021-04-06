function [d_err_est_indiv, d_err_est_comb] = calculate_histogram_distance_err_true_est(err_all, err_est_indiv, err_est_comb, bins, histogram_distance_method)
% Function to calculate histogram distance between true and estimated error
%
% INPUTS:
% err_all: true error
% err_est_indv, err_est_comb: estimated error for individual and combined schemes
% bins: bins for calculating histogram
% histogram_distance_method: type of distance method (e.g. 'total_variation_distance')
%
% OUTPUTS:
% d_err_est_indiv, d_err_est_comb: scaled histogram distance for individual and combined schemes
%
% AUTHOR:
% Lalit Rajendran (lrajendr@purdue.edu)

    % number of individual methods
    num_individual_methods = size(err_est_indiv, 2); 
    % number of combined methods
    num_combined_methods = size(err_est_comb, 2); 

    % bin width
    bin_width = bins(2) - bins(1);
    % individual methods
    for method_index = 1:num_individual_methods
        % d_err_est_indiv(method_index) = calculate_histogram_distance(err_all, err_est_indiv(method_index, :), bins, histogram_distance_method) * bin_width;
        d_err_est_indiv(method_index) = calculate_histogram_distance(err_all, err_est_indiv(:, method_index), bins, histogram_distance_method) * bin_width;
    end

    % combined methods
    for method_index = 1:num_individual_methods
        % d_err_est_comb(method_index) = calculate_histogram_distance(err_all, err_est_comb(method_index, :), bins, histogram_distance_method) * bin_width;
        d_err_est_comb(method_index) = calculate_histogram_distance(err_all, err_est_comb(:, method_index), bins, histogram_distance_method) * bin_width;
    end

end
function [err_est_indiv, err_est_comb] = estimate_error_from_uncertainty_gaussian(err, unc_indiv, unc_comb)
% Function to estimate error distributions from uncertainties
%
% INPUTS:
% err: error
% unc_indiv: individual uncertainties
% unc_comb: combined uncertainties
%
% OUTPUTS:
% err_est_indiv, err_est_comb: estimated errors for indidivual and combined schemes
%
% AUTHOR:
% Lalit Rajendran (lrajendr@purdue.edu)

    % number of individual methods
    num_individual_methods = size(unc_indiv, 2); 
    % number of combined methods
    num_combined_methods = size(unc_comb, 2); 

    % average of the error
    err_avg = mean(err, 'omitnan');

    % individual methods
    for method_index = 1:num_individual_methods
        y = estimate_error_from_uncertainty('gaussian', err_avg, unc_indiv(:, method_index));
        err_est_indiv(:, method_index) = y(:);
    end

    % combined methods
    for method_index = 1:num_combined_methods
        y = estimate_error_from_uncertainty('gaussian', err_avg, unc_comb(:, method_index));
        err_est_comb(:, method_index) = y(:);
    end

end
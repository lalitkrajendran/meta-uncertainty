function [err_rms, unc_indiv_rms, unc_comb_rms, err_avg, confint_err] = calculate_error_uncertainty_rms(err, unc_indiv, unc_comb);
% Function to calculate error and uncertainty statistics such as the mean, rms etc.
%
% INPUTS:
% err: errors
% unc_indiv: individual uncertainties
% unc_comb: combined uncertainties
%
% OUTPUTS:
% err_rms: rms of error
% unc_indiv_rms: rms of individual uncertainty
% unc_comb_rms: rms of combined uncertainty
% err_avg: mean of error
% confint_err: 1/2 of 68% confidence interval 
%
% AUTHOR:
% Lalit Rajendran (lrajendr@purdue.edu)

    % error statistics
    err_rms = rms(err, 'omitnan');
    err_avg = mean(err, 'omitnan');
    % calculate confidence interval for the error
    confint_err = (prctile(err, 84) - prctile(err, 16))/2;
    
    % uncertainty (individual)
    unc_indiv_rms = nanrms(unc_indiv, 1);
    
    % uncertainty (combined)
    unc_comb_rms = nanrms(unc_comb, 1);
end
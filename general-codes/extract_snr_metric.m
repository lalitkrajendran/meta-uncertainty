function snr_metric = extract_snr_metric(SNRmetric, r, c)
% Function to extract the SNR metrics (PPR and MI) for a 
% given point in the FOV
%
% INPUTS:
% SNRmetric: structure from prana results
% r, c: row and column indices of the grid point
%
% OUTPUTS:
% snr_metric: structure containing ppr and mi for this grid point
%
% AUTHORS:
% Lalit Rajendran (lrajendr@purdue.edu)
%
% DATE:
% 2020/06/05

    snr_metric = struct;

    % peak to peak ratio (ppr)
    snr_metric.PPR = SNRmetric.PPR(r, c);

    % mutual information (mi)
    snr_metric.MI = SNRmetric.MI(r, c);
end
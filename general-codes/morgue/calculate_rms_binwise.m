function [err_rms_binwise, sigma_rms_binwise] = calculate_rms_binwise(err_all, sigma_all, max_bin, num_bins)
    
    % generate bins
    edges = linspace(0, max_bin, num_bins);
    % bin the uncertainty
    [~, ~, bin_sigma] = histcounts(sigma_all, edges);

    % initialize binwise rms error and uncertainties
    err_rms_binwise = nans(1, num_bins);
    sigma_rms_binwise = nans(1, num_bins);
    
    % loop through bins and calculate binwise rms
    for bin_index = 1:num_bins
        % find uncertainty values lying in the current bin
        sigma_bin_indices = find(bin_sigma == bin_index);
        % calculate rms of the uncertainties of the values in the
        % current bin
        sigma_rms_binwise(bin_index) = rms(sigma_all(sigma_bin_indices), 'omitnan');
        % calculate rms of the errors of the measurements in the
        % current bin
        err_rms_binwise(bin_index) = rms(err_all(sigma_bin_indices), 'omitnan');
    end
end
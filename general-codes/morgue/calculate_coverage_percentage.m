function coverage = calculate_coverage_percentage(err_all, sigma_all)

    % identify non-nan and valid indices
    non_nan_indices = ~isnan(err_all) & ~isnan(sigma_all);

    % extract non-nan and valid errors and uncertainties
    err_temp = err_all(non_nan_indices);
    if numel(sigma_all) > 1
        sigma_temp = sigma_all(non_nan_indices);
    else
        sigma_temp = sigma_all;
    end

    % calculate coverage (fraction of points with error < uncertainty)
    coverage = sum(abs(err_temp) < abs(sigma_temp))/numel(err_temp) * 100;
end
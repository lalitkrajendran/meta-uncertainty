function weights = calculate_weights_from_resampling(unc_resampling_current, bins)
    weights = cell(1, 3);
    % ===================
    %% unweighted
    % ===================
    combination_method_index = 1;
    
    % calculate weighting matrix
    weights_current = [1, 1, 1];
    weights{combination_method_index}.x = real(weights_current/sum(weights_current));

    weights_current = [1, 1, 1];
    weights{combination_method_index}.y = real(weights_current/sum(weights_current));
    
    % ===================
    %% variance covariance
    % ===================
    combination_method_index = 2;
    
    % calculate covariance matrix
    covariance_matrix_x = cov([unc_resampling_current.imx', unc_resampling_current.mcx', unc_resampling_current.csx']);
    covariance_matrix_y = cov([unc_resampling_current.imy', unc_resampling_current.mcy', unc_resampling_current.csy']);

    % calculate weighting matrix
    weights_current = [sum(1./covariance_matrix_x(:, 1)), sum(1./covariance_matrix_x(:, 2)), sum(1./covariance_matrix_x(:, 3))];
    weights{combination_method_index}.x = real(weights_current/sum(weights_current(:)));

    weights_current = [sum(1./covariance_matrix_y(:, 1)), sum(1./covariance_matrix_y(:, 2)), sum(1./covariance_matrix_y(:, 3))];
    weights{combination_method_index}.y = real(weights_current/sum(weights_current(:)));

    % ===================
    %% entropy    
    % ===================
    combination_method_index = 3;

    entropy = nans(3, 2);
    
    % calculate entropies
    entropy(1, 1) = calculate_shannon_entropy(unc_resampling_current.imx, bins);
    entropy(1, 2) = calculate_shannon_entropy(unc_resampling_current.imy, bins);
    entropy(2, 1) = calculate_shannon_entropy(unc_resampling_current.mcx, bins);
    entropy(2, 2) = calculate_shannon_entropy(unc_resampling_current.mcy, bins);
    entropy(3, 1) = calculate_shannon_entropy(unc_resampling_current.csx, bins);
    entropy(3, 2) = calculate_shannon_entropy(unc_resampling_current.csy, bins);
    
    % calculate weights
    weights_current = 1./entropy(:, 1);
    weights{combination_method_index}.x = real(weights_current/sum(weights_current));

    weights_current = 1./entropy(:, 2);
    weights{combination_method_index}.y = real(weights_current/sum(weights_current));
end
function [weights, pdf_unc_resampling] = calculate_weights_from_resampling_3c(unc_stereo_resampled, bins_x, bins_y, bins_z)
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

    weights_current = [1, 1, 1];
    weights{combination_method_index}.z = real(weights_current/sum(weights_current));
    
    % ===================
    %% variance covariance
    % ===================
    combination_method_index = 2;
    
    % remove outliers
    for method_index = 1:3
        temp = real(unc_stereo_resampled{method_index}.x);
        temp(temp < 0 & temp > bins_x(end)) = 0;
        unc_stereo_resampled{method_index}.x = temp;

        temp = real(unc_stereo_resampled{method_index}.y);
        temp(temp < 0 & temp > bins_y(end)) = 0;
        unc_stereo_resampled{method_index}.y = temp;

        temp = real(unc_stereo_resampled{method_index}.z);
        temp(temp < 0 & temp > bins_z(end)) = 0;
        unc_stereo_resampled{method_index}.z = temp;
    end

    % calculate covariance matrix
    covariance_matrix_x = cov([unc_stereo_resampled{1}.x', unc_stereo_resampled{2}.x', unc_stereo_resampled{3}.x']);
    covariance_matrix_y = cov([unc_stereo_resampled{1}.y', unc_stereo_resampled{2}.y', unc_stereo_resampled{3}.y']);
    covariance_matrix_z = cov([unc_stereo_resampled{1}.z', unc_stereo_resampled{2}.z', unc_stereo_resampled{3}.z']);

    % only retain magnitude of the covariance matrix
    covariance_matrix_x = abs(covariance_matrix_x);
    covariance_matrix_y = abs(covariance_matrix_y);
    covariance_matrix_z = abs(covariance_matrix_z);
    
    % calculate weighting matrix
    weights_current = [sum(1./covariance_matrix_x(:, 1)), sum(1./covariance_matrix_x(:, 2)), sum(1./covariance_matrix_x(:, 3))];
    % % only retain magnitude
    % weights_current = abs(weights_current);
    weights{combination_method_index}.x = real(weights_current/sum(weights_current(:)));

    weights_current = [sum(1./covariance_matrix_y(:, 1)), sum(1./covariance_matrix_y(:, 2)), sum(1./covariance_matrix_y(:, 3))];
    % % only retain magnitude
    % weights_current = abs(weights_current);
    weights{combination_method_index}.y = real(weights_current/sum(weights_current(:)));

    weights_current = [sum(1./covariance_matrix_z(:, 1)), sum(1./covariance_matrix_z(:, 2)), sum(1./covariance_matrix_z(:, 3))];
    % % only retain magnitude
    % weights_current = abs(weights_current);
    weights{combination_method_index}.z = real(weights_current/sum(weights_current(:)));
    
    % ===================
    %% entropy    
    % ===================
    combination_method_index = 3;

    entropy = nans(3, 3);
    pdf_unc_resampling = cell(1, 3);
    % calculate entropies
    for method_index = 1:3
        pdf_unc_resampling{method_index} = struct;
        [entropy(method_index, 1), pdf_unc_resampling{method_index}.x] = calculate_shannon_entropy(unc_stereo_resampled{method_index}.x, bins_x);
        [entropy(method_index, 2), pdf_unc_resampling{method_index}.y] = calculate_shannon_entropy(unc_stereo_resampled{method_index}.y, bins_y);
        [entropy(method_index, 3), pdf_unc_resampling{method_index}.z] = calculate_shannon_entropy(unc_stereo_resampled{method_index}.z, bins_z);        
    end
    
    % calculate weights
    weights_current = 1./entropy(:, 1);
    weights{combination_method_index}.x = real(weights_current'/sum(weights_current));

    weights_current = 1./entropy(:, 2);
    weights{combination_method_index}.y = real(weights_current'/sum(weights_current));

    weights_current = 1./entropy(:, 3);
    weights{combination_method_index}.z = real(weights_current'/sum(weights_current));
end
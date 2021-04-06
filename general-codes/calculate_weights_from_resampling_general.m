function weights = calculate_weights_from_resampling_general(unc_resampled, component_names, bins, combination_methods)
% Function to calculate weights from resampled uncertainties.
%
% INPUTS:
% unc_resampled: resampled uncertainty arrays
% component_names: names of components in the cell array
% bins: histogram bins for entropy calculations
% combination_methods: array of combination/weighting methods
% 
% OUTPUTS:
% weights: array of weights
% pdf_unc_resampling: pdf of resampled uncertainties used in entropy calculations
%
% DEPENDENCIES:
% CompPD package: https://www.jstatsoft.org/article/view/v065i02
% calculate_shannon_entropy: function to calculate shannon entropy
%
% AUTHORS:
% Lalit Rajendran (lrajendr@purdue.edu)
%
% DATE:
% 05/04/2020

    % number of uncertainty methods
    num_uncertainty_methods = numel(unc_resampled);
    % number of resampling trials
    num_resampling_trials = numel(unc_resampled{1}.x);
    % number of components
    num_components = numel(component_names);
    % number of weighting methods
    num_combination_methods = numel(combination_methods);

    % allocate memory for weight matrix
    weights = cell(1, num_combination_methods);

    % ===========================================
    % loop through combination methods and calculate weights
    % ===========================================
    for combination_method_index = 1:num_combination_methods
        % declare structure
        weights{combination_method_index} = struct;
        weights{combination_method_index}.name = combination_methods{combination_method_index};

        % ===================
        %% unweighted
        % ===================            
        if strcmp(combination_methods{combination_method_index}, 'unwt')
            % calculate weighting matrix
            for comp_index = 1:num_components
                weights_current = [1; 1; 1];
                weights{combination_method_index}.(component_names{comp_index}) = real(weights_current/sum(weights_current));
            end
        
        % ===================
        %% variance based weighting
        % ===================            
        elseif contains(combination_methods{combination_method_index}, 'var')
            % loop through components
            for comp_index = 1:num_components
                % allocate memory
                unc_all = nans(num_resampling_trials, num_uncertainty_methods);
                % aggregate
                for uncertainty_method_index = 1:num_uncertainty_methods
                    % aggregate
                    unc_all(:, uncertainty_method_index) = real(unc_resampled{uncertainty_method_index}.(component_names{comp_index}));
                    % remove outliers
                    indices = unc_all(:, uncertainty_method_index) > bins{comp_index}(end);
                    unc_all(indices, uncertainty_method_index) = NaN;
                end            
                                
                % only retain non-nan values
                [r, c] = find(isnan(unc_all));                
                unc_all(r, :) = [];
                if size(unc_all, 1) <= 1
                    fprintf('Empty array\n');
                    weights{combination_method_index}.(component_names{comp_index}) = nans(num_uncertainty_methods, 1);
                    continue;
                end
                
                % % break if there are not enough values to calculate a variance estimate
                % if size(unc_all, 1) <= 100
                %     weights{combination_method_index}.(component_names{comp_index}) = nans(num_uncertainty_methods, 1);
                %     fprintf('insufficient elements to calculate quantiles\n');
                %     continue;
                % end

                % standard covariance calculation
                if strcmp(combination_methods{combination_method_index}, 'var')
                    covariance_matrix = cov(unc_all);
                % projection depth based covariance
                elseif strcmp(combination_methods{combination_method_index}, 'pd-var')
                    try                        
                        % perturb the data by a small amount to ensure it is in 
                        % general position (more than 2 points should not be collinear)
                        unc_all = PertX(unc_all);
                        % estimate number of direction vectors
                        % num_vec = min(size(unc_all, 1)/2, 100);
                        num_vec = size(unc_all, 1) * 10;
                        
                        % obtain the direction vectors for approximate computing the projection depth 
                        % and its associated estimators
                        % u = AppVecAPD(unc_all, num_vec*10);
                        u = ExVecAPDHD(unc_all);
                        
                        % extract interquartile contour
                        vpmat = APC3D(unc_all, u, 0.5, true, false);
                        % vpmat = PC3D(unc_all, u, 0.5, true, false);
                    catch
                        fprintf('interquartile calculation failed\n');
                        weights{combination_method_index}.(component_names{comp_index}) = nans(num_uncertainty_methods, 1);
                        continue;
                    end
                    
                    % calculate covariance of the interquartile contour
                    if numel(vpmat{1}) > 0                        
                        % covariance_matrix = cov([vpmat{1}(:, 1), vpmat{1}(:, 2), vpmat{1}(:, 3)]);
                        covariance_matrix = cov(vpmat{1});
                    else
                        covariance_matrix = nans(num_uncertainty_methods);
                    end
                end

                covariance_matrix = diag(diag(covariance_matrix));
                % calculate inverse of the covariance matrix
                % inv_cov = abs(1./covariance_matrix);
                inv_cov = inv(covariance_matrix);

                % calculate weights
                % weights_current = sum(1./covariance_matrix, 1);
                % weights_current = inv_cov * ones(num_uncertainty_methods, 1) ./ sum(inv_cov(:));
                weights_current = sum(inv_cov, 2) / sum(inv_cov, 'all');
                weights_current = abs(weights_current);
                weights_current = weights_current/sum(weights_current);                

                % normalize and assign weights
                weights{combination_method_index}.(component_names{comp_index}) = real(weights_current/sum(weights_current(:)));
            end
        
        % ===================
        %% entropy based weighting
        % ===================
        elseif strcmp(combination_methods{combination_method_index}, 'entropy')
            
            % allocate memory
            entropy = nans(num_uncertainty_methods, num_components);
            pdf_unc_resampling = cell(1, num_uncertainty_methods);            
            % loop through components
            for comp_index = 1:num_components
                % allocate memory
                information = nans(num_uncertainty_methods, num_uncertainty_methods);
                % calculate entropies
                for uncertainty_method_index = 1:num_uncertainty_methods
                    pdf_unc_resampling{uncertainty_method_index} = struct;
                    % calculate entropies
                    [entropy(uncertainty_method_index, comp_index), pdf_unc_resampling{uncertainty_method_index}.(component_names{comp_index})] = ...
                                            calculate_shannon_entropy(unc_resampled{uncertainty_method_index}.(component_names{comp_index}), bins{comp_index});                    

                    % % calculate individual entropies
                    % information(uncertainty_method_index, uncertainty_method_index) = calculate_self_entropy(unc_resampled{uncertainty_method_index}.(component_names{comp_index}), bins{comp_index});
                    % % calculate mutual entropies
                    % for method_index = 1:num_uncertainty_methods
                    %     if method_index == uncertainty_method_index
                    %         continue;
                    %     else
                    %         unc1 = unc_resampled{uncertainty_method_index}.(component_names{comp_index});
                    %         unc2 = unc_resampled{method_index}.(component_names{comp_index});

                    %         information(uncertainty_method_index, method_index) = calculate_mutual_information(unc1, unc2, bins, bins);
                    %     end
                    % end
                end

                % calculate weights
                weights_current = 1./entropy(:, comp_index);
                
                % % calculate inverse of the entropy matrix
                % inv_inf = inv(information);
                
                % % calculate weights from the inverse
                % weights_current = sum(inv_inf, 2)/sum(inv_inf, 'all');
                % weights_current = abs(weights_current);

                % normalize and assign weights
                weights{combination_method_index}.(component_names{comp_index}) = real(weights_current/sum(weights_current(:)));
            end
        end
    end
end
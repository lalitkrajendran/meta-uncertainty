function [valid_elements_all, err_valid, unc_indiv_valid, unc_comb_valid] = extract_valid_errors_uncertainties(err, unc_indiv, unc_comb, components, min_error_threshold, max_error_threshold)
% Function to extract valid errors and uncertainties
%
% INPUTS:
% err: errors
% unc_indiv: individual uncertainties
% unc_comb: combined uncertainties
% components: displacement components ('x', 'y' etc.)
%
% OUTPUTS:
% valid_elements_all: logical array of valid elements (1 for valid, 0 for invalid)
% err_valid: valid errors
% unc_indiv_valid: valid individual uncertainties
% unc_comb_valid: valid combined uncertainties
%
% AUTHOR:
% Lalit Rajendran (lrajendr@purdue.edu)

    % number of components
    num_components = numel(components);
    % number of individual uncertainty methods
    num_individual_methods = size(unc_indiv.(components{1}), 2);
    % number of combined uncertainty methods
    num_combined_methods = size(unc_comb.(components{1}), 2);

    % ==========================
    % identify valid measurements
    % ==========================
    for component_index = 1:num_components
        % component name
        component_name = components{component_index};
        % error            
        valid_elements(:, 1) = identify_valid_elements(abs(err.(component_name)), min_error_threshold(component_index), max_error_threshold(component_index));

        % uncertainty (individual)
        for method_index = 1:num_individual_methods
            valid_elements(:, 1 + method_index) = identify_valid_elements(abs(unc_indiv.(component_name)(:, method_index)), min_error_threshold(component_index), max_error_threshold(component_index));
        end

        % combined (individual)
        for method_index = 1:num_combined_methods
            valid_elements(:, 1 + num_individual_methods + method_index) = identify_valid_elements(abs(unc_comb.(component_name)(:, method_index)), min_error_threshold(component_index), max_error_threshold(component_index));
        end

        % find valid elements across all errors and uncertainties for this component
        valid_elements_component(:, component_index) = prod(valid_elements, 2);
    end

    % find elements valid across both components
    valid_elements_all = logical(prod(valid_elements_component, 2));
    num_valid_elements = sum(double(valid_elements_all));

    % ==========================
    % extract valid measurements
    % ==========================    
    for component_index = 1:num_components
        % component name
        component_name = components{component_index};

        % error
        err_valid.(component_name) = err.(component_name)(valid_elements_all');

        % uncertainty (individual)
        for method_index = 1:num_individual_methods
            unc_indiv_valid.(component_name)(:, method_index) = unc_indiv.(component_name)(valid_elements_all, method_index);
        end

        % uncertainty (combined)
        for method_index = 1:num_combined_methods
            unc_comb_valid.(component_name)(:, method_index) = unc_comb.(component_name)(valid_elements_all, method_index);
        end
    end
end
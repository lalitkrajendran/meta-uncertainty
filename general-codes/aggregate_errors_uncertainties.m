function [err_all, unc_indiv_all, unc_comb_all] = aggregate_errors_uncertainties(err, unc_indiv, unc_comb, components)
% Function to aggregrate errors and uncertainties across components
%
% INPUTS:
% err: errors
% unc_indiv: individual uncertainties
% unc_comb: combined uncertainties
% components: displacement components ('x', 'y' etc.)
% 
% OUTPUTS:
% err_all: aggregated errors
% unc_indiv_all: aggregated individual uncertainties
% unc_comb_all: aggregated combined uncertainties
%
% AUTHOR:
% Lalit Rajendran (lrajendr@purdue.edu)

    % number of components
    num_components = numel(components);
    % number of elements in each component
    num_elements = numel(err.(components{1}));
    % number of individual uncertainty methods
    num_individual_methods = size(unc_indiv.(components{1}), 2);
    % number of combined uncertainty methods
    num_combined_methods = size(unc_comb.(components{1}), 2);
    
    % ==========================
    % aggregate valid measurements
    % ==========================
    err_all = nans(num_elements * num_components, 1);
    unc_indiv_all = nans(num_elements * num_components, num_individual_methods);    
    unc_comb_all = nans(num_elements * num_components, num_combined_methods);    
    
    % loop through components
    for component_index = 1:num_components
        % component name
        component_name = components{component_index};
        % start index
        start_index = (component_index - 1) * num_elements + 1;
        stop_index = component_index * num_elements;

        % error
        err_all(start_index:stop_index, 1) = err.(component_name);

        % uncertainty (individual)
        for method_index = 1:num_individual_methods
            unc_indiv_all(start_index:stop_index, method_index) = unc_indiv.(component_name)(:, method_index);
        end

        % uncertainty (combined)
        for method_index = 1:num_combined_methods
            unc_comb_all(start_index:stop_index, method_index) = unc_comb.(component_name)(:, method_index);
        end
    end

end


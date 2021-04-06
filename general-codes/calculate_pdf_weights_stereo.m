function pdf_w = calculate_pdf_weights_stereo(wt, bins_w, valid_trials, num_individual_methods, num_combined_methods, components)
% Function to calculate pdf of weights from resampling
%
% INPUTS:
% wt: cell array containing weights
% bins_w: bins to use for calculating pdf
% valid_trials: logical array of valid elements
% num_individual_methods: number of individual methods
% num_combined_methods: number of combined methods
% components: cell array of displacement components ('x', 'y', etc.)
%
% OUTPUTS:
% pdf_w: pdf of weights: cell array of same size as wt
%
% AUTHOR:
% Lalit Rajendran (lrajendr@purdue.edu)
        
    % number of components
    num_components = numel(components);
    pdf_w = struct;

    % create cell array to hold results
    for component_index = 1:num_components
        component_name = components{component_index};
        pdf_w.(component_name) = cell(num_combined_methods, num_individual_methods);
        % loop through combined methods
        for combined_method_index = 1:num_combined_methods
            % loop through individual methods
            for individual_method_index = 1:num_individual_methods
                % weights
                wt_temp = wt{combined_method_index}.(component_name)(valid_trials, individual_method_index);
                % calculate pdf        
                pdf_w.(component_name){combined_method_index, individual_method_index} = histcounts(wt_temp(isfinite(wt_temp)), bins_w, 'normalization', 'pdf');
            end
        end    
    end
end
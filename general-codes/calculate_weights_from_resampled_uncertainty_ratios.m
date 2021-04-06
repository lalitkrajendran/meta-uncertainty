function [wt, unc_ratio_metric, unc_ratio_fit] = calculate_weights_from_resampled_uncertainty_ratios(unc_indiv, unc_resampled, metric_name, individual_method_array, ...
                                    percentage_particles_remove_array, components, num_resampling_trials, data_type)
% Function to calculate weights for uncertainty schemes based on the variation with particle addition/removal
%
% INPUTS:
% unc_indiv: individual uncertainties
% unc_comb: combined uncertainties
% metric_name: statistical metric to condense resampling results (e.g. 'iqr')
% individual_method_array: names of individual uncertainty methods
% percentage_particles_remove_array: array of particle addition/removal
% component: co-ordinates ('x', 'y')
% num_resampling_trials: number of resampling trials
% data_type: data type for storing the resampled uncertainties
%
% OUTPUTS:
% wt: weights
% unc_ratio_metric: statistics of the uncertainty ratios
% unc_ratio_fit: straight line fit to the uncertainty ratio statistics vs particle perturbation
%
% AUTHOR:
% Lalit Rajendran (lrajendr@purdue.edu)

    % number of individual methods
    num_individual_methods = numel(individual_method_array);
    % number of particle perturbation cases
    num_ppr = numel(percentage_particles_remove_array);
    % number of components
    num_components = numel(components);
    
    % ==========================
    % initialize
    % ==========================
    unc_ratio_metric = struct;
    unc_ratio_fit = struct;
    wt = struct;
    
    % ==========================
    % loop through components
    % ==========================
    for component_index = 1:num_components
        % component name
        component_name = components{component_index};
        % ==========================
        % loop through uncertainty methods
        % ==========================
        for individual_method_index = 1:num_individual_methods
            % method name
            method_name = lower(individual_method_array{individual_method_index});

            % extract individual uncertainty
            unc_sub = unc_indiv{individual_method_index}.(component_name);

            % initialze arrays
            unc_ratio = nans(num_resampling_trials, num_ppr);

            for particle_remove_index = 1:num_ppr
                % current percentage to remove
                percentage_particles_remove = percentage_particles_remove_array(particle_remove_index);

                if strcmp(data_type, 'struct')
                    unc_ratio(:, particle_remove_index) = unc_resampled{particle_remove_index}.([method_name component_name]) / unc_sub;
                elseif strcmp(data_type, 'cell')
                    unc_ratio(:, particle_remove_index) = unc_resampled{particle_remove_index}{individual_method_index}.(component_name) / unc_sub;
                end
            end

            % --------------------------
            % calculate iqr of the ratio across the resampling trials
            % --------------------------
            if strcmpi(metric_name, 'iqr')
                unc_ratio_metric.(component_name)(individual_method_index, :) = iqr(unc_ratio, 1);
            elseif strcmpi(metric_name, 'rms')
                unc_ratio_metric.(component_name)(individual_method_index, :) = nanrms(unc_ratio, 1);
            elseif strcmpi(metric_name, 'std')
                unc_ratio_metric.(component_name)(individual_method_index, :) = std(unc_ratio, 1, 'omitnan');
            end

            % --------------------------
            % calculate rate of change
            % --------------------------
            % extract values
            f = unc_ratio_metric.(component_name)(individual_method_index, :)';

            % straight line fit
            fitobject = fit(percentage_particles_remove_array', squeeze(f), 'poly1');
            
            % extract fit values
            unc_ratio_fit.(component_name)(individual_method_index, :) = fitobject(percentage_particles_remove_array);

            % calculate weights
            wt.(component_name)(individual_method_index) = 1 ./ abs(fitobject.p1);
        end

        % normalize weights so they add up to 1
        wt.(component_name) = wt.(component_name) / sum(wt.(component_name));
    end
end
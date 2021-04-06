function plot_unc_ratio_with_particle_removal(unc_ratio, unc_ratio_fit, trial_index, ppr, uncertainty_method_names, resampling_method_names, metric_name, ...
                                                user_screen_resolution)

    % number of methods
    num_uncertainty_methods = numel(uncertainty_method_names);
    num_resampling_methods = numel(resampling_method_names);

    colors = lines(num_uncertainty_methods);
    
    figure
    for resampling_method_index = 1:num_resampling_methods
        subplot(1, num_resampling_methods, resampling_method_index)
        for uncertainty_method_index = 1:num_uncertainty_methods
            plot(ppr, squeeze(unc_ratio{resampling_method_index}.x(trial_index, uncertainty_method_index, :)), '*', 'color', colors(uncertainty_method_index, :))
            hold on
            h = plot(ppr, squeeze(unc_ratio_fit{resampling_method_index}.x(trial_index, uncertainty_method_index, :)), 'color', [colors(uncertainty_method_index, :), 0.5]);
            if resampling_method_index == 1
                l(uncertainty_method_index) = h;
            end
            
            plot(ppr, squeeze(unc_ratio{resampling_method_index}.y(trial_index, uncertainty_method_index, :)), 'o', 'color', colors(uncertainty_method_index, :))
            plot(ppr, squeeze(unc_ratio_fit{resampling_method_index}.y(trial_index, uncertainty_method_index, :)), '--', 'color', [colors(uncertainty_method_index, :), 0.5]);

        end    
        
        box off
        ylim([0, 3])

        % annotate figure
        xlabel('Particle Perturbation %')
        if resampling_method_index == 1
            ylabel(metric_name)
            % plot invisible markers for x and y
            subplot(1, num_resampling_methods, 1)
            l_x = plot(NaN, NaN, 'k*');
            l_y = plot(NaN, NaN, 'ko');
            
            % add legend
            lgd = legend([l, l_x, l_y], {uncertainty_method_names{:}, 'x', 'y'}, 'location', 'northoutside', 'orientation', 'horizontal');                        
            
            % set figure position        
            set(gcf, 'resize', 'off');
            set(gcf, 'units', 'inches', 'position', [267   377   972   364]/user_screen_resolution);
            set(gcf, 'resize', 'off');
            lgd.Position(1:2) = [0.3, 0.9];
        end
        title(resampling_method_names{resampling_method_index})            

        pos = get(gca, 'position');
        pos(2) = 0.15;
        pos(4) = 0.65;
        set(gca, 'position', pos)                    
    end

    drawnow();
end
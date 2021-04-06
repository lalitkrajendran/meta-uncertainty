function plot_unc_ratio_with_particle_removal_stereo(unc_ratio, unc_ratio_fit, ppr, uncertainty_method_names, resampling_method_names, metric_name, ...
                                                    component_names, line_symbols, symbols, user_screen_resolution)

    % number of methods
    num_uncertainty_methods = numel(uncertainty_method_names);
    num_resampling_methods = numel(resampling_method_names);
    num_components = numel(component_names);
    colors = lines(1); %num_uncertainty_methods);
    component_names_plot = {'U'; 'V'; 'W'};
    
    figure
    for component_index = 1:num_components
        for resampling_method_index = 1:num_resampling_methods
            % subplot_index = (component_index - 1) * num_resampling_methods + resampling_method_index;
            subplot_index = (resampling_method_index - 1) * num_components + component_index;
            subplot(num_resampling_methods, num_components, subplot_index)
            for uncertainty_method_index = 1:num_uncertainty_methods
                % plot(ppr, squeeze(unc_ratio{resampling_method_index}.(component_names{component_index})(uncertainty_method_index, :)), ...
                %                     '*', 'color', colors(uncertainty_method_index, :))
                plot(ppr, squeeze(unc_ratio{resampling_method_index}.(component_names{component_index})(uncertainty_method_index, :)), ...
                                    symbols{uncertainty_method_index}, 'color', colors(1, :))
                hold on
                % h = plot(ppr, squeeze(unc_ratio_fit{resampling_method_index}.(component_names{component_index})(uncertainty_method_index, :)), ...
                %                     'color', [colors(uncertainty_method_index, :), 0.5]);
                h = plot(ppr, squeeze(unc_ratio_fit{resampling_method_index}.(component_names{component_index})(uncertainty_method_index, :)), ...
                                    line_symbols{uncertainty_method_index}, 'color', [colors(1, :), 0.5]);
                if subplot_index == 1 %resampling_method_index == 1
                    l(uncertainty_method_index) = h;
                end
            end    
    
            box off
            % ylim([0, 3])
            ylim([0 1])
    
            % annotate figure
            if resampling_method_index == 1
                % title(resampling_method_names{resampling_method_index})            
                title(component_names_plot{component_index})
            end

            if resampling_method_index == num_resampling_methods
                xlabel('Particle Perturbation %')
            else
                set(gca, 'xticklabel', []);
            end
            if component_index == 1
                ylabel([metric_name]) % ', \sigma_' component_names_plot{component_index}])
            else
                set(gca, 'yticklabel', [])
            end
            
            
            pos = get(gca, 'position');
            pos(2) = 0.15;
            pos(4) = 0.65;
            set(gca, 'position', pos)                    
        end
    end


    % add legend
    % lgd = legend([l, l_x, l_y], {uncertainty_method_names{:}, 'x', 'y'}, 'location', 'northoutside', 'orientation', 'horizontal');                        
    lgd = legend(l, {uncertainty_method_names{:}}, 'location', 'northoutside', 'orientation', 'horizontal', 'fontsize', 14);                        

    % set figure position        
    set(gcf, 'resize', 'off');
    % set(gcf, 'units', 'inches', 'position', [267   377   972   900]/user_screen_resolution);
    set(gcf, 'units', 'inches', 'position', [254   501   850   350]/user_screen_resolution);
    
    set(gcf, 'resize', 'off');
    lgd.Position(1:2) = [0.35, 0.9];


    drawnow();
end
function lgd = plot_unc_ratio_with_particle_removal_single(unc_ratio, unc_ratio_fit, weights, trial_index, ppr, uncertainty_method_names, resampling_method_names, metric_name, ...
    symbols, line_symbols, user_screen_resolution)

    % number of methods
    num_uncertainty_methods = numel(uncertainty_method_names);
    num_resampling_methods = numel(resampling_method_names);

    % colors = lines(num_uncertainty_methods);
    colors = lines(1);

    component_names = {'x'; 'y'};
    num_components = numel(component_names);

    figure
    for component_index = 1:num_components
        component_name = component_names{component_index};
        subplot(1, num_components, component_index)
        for uncertainty_method_index = 1:num_uncertainty_methods
            % plot(ppr, squeeze(unc_ratio.(component_name)(trial_index, uncertainty_method_index, :)), '*', 'color', colors(uncertainty_method_index, :))
            % hold on
            % h = plot(ppr, squeeze(unc_ratio_fit.(component_name)(trial_index, uncertainty_method_index, :)), 'color', [colors(uncertainty_method_index, :), 0.5]);
            
            h1 = plot(ppr, squeeze(unc_ratio.(component_name)(trial_index, uncertainty_method_index, :)), symbols{uncertainty_method_index}, 'color', colors(1, :));
            l1(uncertainty_method_index) = h1;
            hold on
            h = plot(ppr, squeeze(unc_ratio_fit.(component_name)(trial_index, uncertainty_method_index, :)), line_symbols{uncertainty_method_index}, 'color', [colors(1, :), 0.5]);
            l(uncertainty_method_index) = h;
            
            % add text to indicate weights
            tx = ppr(end) - 2.5;
            ty = unc_ratio_fit.(component_name)(trial_index, uncertainty_method_index, end) + 0.05;
            str = num2str(weights.(component_name)(trial_index, uncertainty_method_index), '%.2f');
            % text(tx, ty, str, 'color', colors(uncertainty_method_index, :), 'fontweight', 'bold')
            text(tx, ty, str, 'color', colors(1, :), 'fontweight', 'bold')
            
            % plot(ppr, squeeze(unc_ratio.y(trial_index, uncertainty_method_index, :)), 'o', 'color', colors(uncertainty_method_index, :))
            % plot(ppr, squeeze(unc_ratio_fit.y(trial_index, uncertainty_method_index, :)), '--', 'color', [colors(uncertainty_method_index, :), 0.5]);
        end    
    
        box off
        ylim([0 0.75])

        % annotate figure
        xlabel('Particle Addition %')
        ylabel([metric_name ', ' upper(component_name)])
        % if component_index == 1
        %     ylabel(metric_name)
        % else
        %     set(gca, 'yticklabels', []);
        % end
        
        % title(upper(component_name))

    end

    % add legend
    % lgd = legend(l, {uncertainty_method_names{:}}, 'location', 'northoutside', 'orientation', 'horizontal');                        
    lgd = legend(l1, {uncertainty_method_names{:}}, 'location', 'northoutside', 'orientation', 'horizontal');                        

    
    set(gcf, 'resize', 'off');
    set(gcf, 'units', 'inches', 'position', [388   489   665   315]/user_screen_resolution)
    set(gcf, 'resize', 'off');

    % set_subplots_height(gcf, 0.65);
    set_subplots_position_y(gcf, 0.15, 0.65);
    lgd.Position(1:2) = [0.35, 0.9];
    % if resampling_method_index == 1
    %     ylabel(metric_name)
    %     % plot invisible markers for x and y
    %     subplot(1, num_resampling_methods, 1)
    %     l_x = plot(NaN, NaN, 'k*');
    %     l_y = plot(NaN, NaN, 'ko');

    %     % add legend
    %     lgd = legend([l, l_x, l_y], {uncertainty_method_names{:}, 'x', 'y'}, 'location', 'northoutside', 'orientation', 'horizontal');                        

    %     % set figure position        
    %     set(gcf, 'resize', 'off');
    %     set(gcf, 'units', 'inches', 'position', [267   377   972   364]/user_screen_resolution);
    %     set(gcf, 'resize', 'off');
    %     lgd.Position(1:2) = [0.3, 0.9];
    %     end
    % title(resampling_method_names{resampling_method_index})            

    % pos = get(gca, 'position');
    % pos(2) = 0.15;
    % pos(4) = 0.65;
    % set(gca, 'position', pos)                    
    

    drawnow();
end
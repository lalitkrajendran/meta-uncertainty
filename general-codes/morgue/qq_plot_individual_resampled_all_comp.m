function qq_plot_individual_resampled_all_comp(sigma_valid_individual, sigma_valid_resampled, ...
    individual_method_array, max_error_threshold, component_names, num_skip, colors, symbols, user_screen_resolution)

    num_components = size(sigma_valid_individual, 2);
    num_individual_methods = numel(individual_method_array);

    % plot marker size
    marker_size = 4;

    figure
    % adjust figure position
    set(gcf, 'resize', 'off');
    drawnow();
    % set(gcf, 'units', 'inches', 'Position', [137         545        1402         355]/user_screen_resolution)
    % set(gcf, 'units', 'inches', 'Position', [137   589   908   311]/user_screen_resolution)
    set(gcf, 'units', 'inches', 'Position', [50   50   908   900]/user_screen_resolution)
    drawnow();
    
    for comp_index = 1:num_components
        for method_index = 1:num_individual_methods
            subplot_index = (comp_index - 1) * num_individual_methods + method_index;
            subplot(num_components, num_individual_methods, subplot_index)
            % plot 1:1 line
            plot([0 max_error_threshold(comp_index)], [0 max_error_threshold(comp_index)], 'k');
            hold on
            
            % extract plot data for this method
            % method name
            method_name = individual_method_array{method_index};
            % individual uncertainty
            x = abs(sigma_valid_individual{method_index, comp_index})';
            % resampled uncertainty
            y = abs(sigma_valid_resampled{method_index, comp_index})';
            % marker color
            color_current = colors(1, :);
            
            % make plot
            l = qqplot(x(1:num_skip:end), y(1:num_skip:end));
            
            % adjust plot
            box off
            axis equal    
            % axis([0 0.1 0 0.1])
            axis([0 max_error_threshold(comp_index) 0 max_error_threshold(comp_index)])

            % adjust line properties
            set(l(1), 'marker', 'o');
            set(l(1), 'markeredgecolor', color_current);
            set(l(1), 'linewidth', 0.2);
            set(l(1), 'markersize', marker_size);
            set(l(2), 'visible', 'off');
            set(l(3), 'visible', 'off');
                        
            % annotate figure
            if method_index == 1 && comp_index == 1
                xlabel('Individual (pix.)'); %Uncertainty (pix.)'); 
                ylabel('Resampled (pix.)'); %Uncertainty (pix.)'); 
            else
                xlabel('');
                ylabel('');

                set(gca, 'yticklabel', []);
                set(gca, 'xticklabel', []);
            end
            if comp_index == 1
                title(method_name);
            end
            if method_index == 1
                text(-max_error_threshold(comp_index)*0.75, max_error_threshold(comp_index)/2, component_names{comp_index}, 'fontsize', 14, 'fontweight', 'bold')                
            end
        end
    end

    
    % adjust figure position
    set(gcf, 'resize', 'off');
    drawnow();
    % set(gcf, 'units', 'inches', 'Position', [137         545        1402         355]/user_screen_resolution)
    % set(gcf, 'units', 'inches', 'Position', [137   589   908   311]/user_screen_resolution)
    set(gcf, 'units', 'inches', 'Position', [50   50   908   900]/user_screen_resolution)
    drawnow();
end
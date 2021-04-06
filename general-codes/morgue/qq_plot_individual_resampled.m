function qq_plot_individual_resampled(sigma_valid_individual, sigma_valid_resampled, ...
    individual_method_array, max_error_threshold, num_skip, colors, symbols, user_screen_resolution)

    num_individual_methods = numel(individual_method_array);

    % plot marker size
    marker_size = 4;

    % -----------------------
    % individual methods
    % -----------------------
    qq_indiv = [];
    for method_index = 1:num_individual_methods
        subplot(1, num_individual_methods, method_index)
        % plot 1:1 line
        plot([0 max_error_threshold], [0 max_error_threshold], 'k');
        hold on
        
        % extract plot data for this method
        % method name
        method_name = individual_method_array{method_index};
        % individual uncertainty
        x = abs(sigma_valid_individual{method_index})';
        % resampled uncertainty
        y = abs(sigma_valid_resampled{method_index})';
        % marker color
        color_current = colors(1, :);
        
        % make plot
        l = qqplot(x(1:num_skip:end), y(1:num_skip:end));
        
        % adjust plot
        box off
        axis equal    
        % axis([0 0.1 0 0.1])
        axis([0 max_error_threshold 0 max_error_threshold])

        % adjust line properties
        set(l(1), 'marker', 'o');
        set(l(1), 'markeredgecolor', color_current);
        set(l(1), 'linewidth', 0.2);
        set(l(1), 'markersize', marker_size);
        set(l(2), 'visible', 'off');
        set(l(3), 'visible', 'off');
        
        % collect line handle
        qq_indiv = [qq_indiv, l(1)];
        
        % annotate figure
        if method_index == 1
            xlabel('Individual Uncertainty (pix.)'); 
            ylabel('Resampled Uncertainty (pix.)'); 
        else
            xlabel('');
            ylabel('');

            set(gca, 'yticklabel', []);
            set(gca, 'xticklabel', []);
        end
        title(method_name);
    end

    
    % adjust figure position
    set(gcf, 'resize', 'off');
    drawnow();
    % set(gcf, 'units', 'inches', 'Position', [137         545        1402         355]/user_screen_resolution)
    set(gcf, 'units', 'inches', 'Position', [137   589   908   311]/user_screen_resolution)
    drawnow();
end
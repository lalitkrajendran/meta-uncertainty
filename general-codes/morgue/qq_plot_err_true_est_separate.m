function qq_plot_err_true_est_separate(err_all_valid, err_est_valid_individual, err_est_valid_combined, ...
    individual_method_array, combination_method_array, max_error_threshold, num_skip, colors, symbols, user_screen_resolution)

    num_individual_methods = numel(individual_method_array);
    num_combination_methods = numel(combination_method_array);

    % extract error values
    x = abs(err_all_valid)';
    x = x(isfinite(x));

    % plot marker size
    marker_size = 4;

    % -----------------------
    % individual methods
    % -----------------------
    qq_indiv = [];
    for method_index = 1:num_individual_methods+num_combination_methods
        subplot(1, num_individual_methods + num_combination_methods, method_index)
        % plot 1:1 line
        plot([0 max_error_threshold], [0 max_error_threshold], 'k');
        hold on
        
        % extract plot data for this method
        if method_index <= num_individual_methods
            % method name
            method_name = individual_method_array{method_index};
            % estimated error
            y = abs(err_est_valid_individual{method_index})';
            % marker color
            color_current = colors(1, :);
        else
            % method name
            method_name = combination_method_array{method_index-3};
            % estimated error
            y = abs(err_est_valid_combined{method_index-3})';
            % marker color
            color_current = colors(2, :);
        end
        
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
            xlabel('True Error (pix.)'); 
            ylabel('Estimated Error (pix.)'); 
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
    set(gcf, 'units', 'inches', 'Position', [137         545        1402         355]/user_screen_resolution)
    drawnow();
end
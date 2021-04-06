function lgd = qq_plot_err_true_est_separate_all_comp_02(err_all_valid, err_est_valid_individual, err_est_valid_combined, ...
    individual_method_array, combination_method_array, component_names, max_error_threshold, num_skip, colors, symbols, user_screen_resolution)

    num_components = size(err_all_valid, 2);
    num_individual_methods = numel(individual_method_array);
    num_combined_methods = numel(combination_method_array);
    num_methods_total = num_individual_methods+num_combined_methods;
    % plot marker size
    marker_size = 4;    

    % -----------------------
    % loop through methods
    % -----------------------
    qq_rms_indiv = nans(num_components, num_individual_methods);
    qq_rms_comb = nans(num_components, num_combined_methods);
    for component_index = 1:num_components
        % extract error values
        x = err_all_valid{component_index}';
        x = x(isfinite(x));        
        l_all = [];
        subplot(1, num_components, component_index)
        plot([-max_error_threshold(component_index), max_error_threshold(component_index)], [-max_error_threshold(component_index), max_error_threshold(component_index)], 'k')
        hold on
        % loop through method types (individual and combined)
        for method_type = 1:2
            % extract individual method data
            if method_type == 1
                % extract errors
                err_est = err_est_valid_individual;
                % number of methods
                num_methods = num_individual_methods;
            % extract combined method data
            else
                % extract errors
                err_est = err_est_valid_combined;
                % number of methods
                num_methods = num_combined_methods;
            end

            % loop through method names
            for method_index = 1:num_methods
                % extract estimated error
                % y = err_est{method_index}';
                if num_methods > 1
                    y = err_est{method_index, component_index};
                else
                    y = err_est{component_index};
                end
                
                % plot
                l = qqplot(x(1:num_skip:end), y(1:num_skip:end));
                % l = qqplot(abs(x(1:num_skip:end)), abs(y(1:num_skip:end)));
                hold on
                % calculate rms deviation
                if method_type == 1
                    qq_rms_indiv(component_index, method_index) = rms(l(1).YData - l(1).XData, 'omitnan');
                else                
                    qq_rms_comb(component_index, method_index) = rms(l(1).YData - l(1).XData, 'omitnan');
                end

                % adjust plot
                set(l(1), 'marker', symbols{method_index});
                set(l(1), 'markeredgecolor', colors(method_type, :));
                set(l(1), 'linewidth', 0.2);
                set(l(1), 'markersize', marker_size);
                
                set(l(2), 'visible', 'off');
                set(l(3), 'visible', 'off');
                
                % aggregate handles to lines
                l_all = [l_all, l(1)];                        
            end
        end

        % adjust plot
        box off
        axis equal
        if component_index < num_components - 1
            set(gca, 'xtick', -0.4:0.2:0.4)
            set(gca, 'ytick', -0.4:0.2:0.4)
        end    
        set(gca, 'YAxisLocation', 'origin')
        set(gca, 'XAxisLocation', 'origin')
        axis([-max_error_threshold(component_index), max_error_threshold(component_index), -max_error_threshold(component_index), max_error_threshold(component_index)])
        xl = xlabel('\epsilon_{True} (pix.)', 'fontsize', 10); %16)
        set(xl, 'HorizontalAlignment', 'center');
        yl = ylabel('\epsilon_{Est} (pix.)', 'fontsize', 10); %16)
        set(yl, 'HorizontalAlignment', 'right');
        title(component_names{component_index});
        drawnow();
    end

    % l_all
    % {individual_method_array{:}, combination_method_array{:}}
    lgd = legend(l_all, {individual_method_array{:}, combination_method_array{:}}, 'location', 'northoutside', 'orientation', 'horizontal'); %, 'NumColumns', 2);
    lgd.Position = [0.4000 0.8096 0.2077 0.0412];
    % adjust figure size
    set(gcf, 'resize', 'off');
    drawnow();
    % set(gcf, 'units', 'inches', 'Position', [500   423   1000   518]/user_screen_resolution);
    set(gcf, 'units', 'inches', 'Position', [23   311   1656   600]/user_screen_resolution);
    drawnow();

end
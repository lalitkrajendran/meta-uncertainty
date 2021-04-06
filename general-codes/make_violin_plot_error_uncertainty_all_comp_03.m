function violins = make_violin_plot_error_uncertainty_all_comp_03(err, sigma_indv, sigma_comb, err_rms, sigma_rms_indv, sigma_rms_comb, bins, ...
                                                                    individual_method_array, combination_method_array, combination_method_plot_index_array, ...
                                                                    component_names, component_names_plot, colors, user_screen_resolution, max_error_threshold)
    
    % --------------------------
    % calculate number of variables involved
    % --------------------------
    num_components = numel(component_names);
    num_individual_methods = numel(individual_method_array);
    num_combined_methods = numel(combination_method_plot_index_array);
    num_methods_total = num_individual_methods + num_combined_methods + 1;

    figure
    drawnow();    
    
    % --------------------------
    % loop through components
    % --------------------------
    for component_index = 1:num_components
        % extract component name
        component_name = component_names{component_index};
        
        subplot(num_components, 1, component_index)
        % calculate number of error values
        N_max = numel(err.(component_name)); 

        % --------------------------
        % make violin plot arrays and colors
        % --------------------------
        plot_array = nans(N_max, num_individual_methods + num_combined_methods + 1);
        colors_array = nans(num_methods_total, 3);
        % individual methods
        for method_index = 1:num_individual_methods
            plot_array(:, method_index) = sigma_indv.(component_name)(:, method_index);
            colors_array(method_index, :) = colors(1, :);
        end
        % error
        plot_array(:, num_individual_methods+1) = abs(err.(component_name));
        colors_array(num_individual_methods+1, :) = [0, 0, 0];
        % combined method
        for method_index = 1:num_combined_methods
            plot_array(:, method_index + num_individual_methods + 1) = sigma_comb.(component_name)(:, combination_method_plot_index_array(method_index));
            colors_array(method_index + num_individual_methods + 1, :) = colors(2, :);
        end

        % plot violins
        violins = violin_plot_pdf(plot_array, bins{component_index}, {individual_method_array{:}, 'Error', combination_method_array{combination_method_plot_index_array}}, colors_array);

        % face colors
        color_all = cell(1, numel(violins));

        % --------------------------
        % adjust violin colors and add rms info
        % --------------------------        
        for violin_index = 1:numel(violins)
            % x co-ordinates of current violin plot
            x = violins(violin_index).XData;

            % individual uncertainty methods
            if violin_index <= num_individual_methods
                % set color
                violins(violin_index).FaceColor = colors(1, :);                
                % extract rms
                y = sigma_rms_indv(violin_index, component_index);
            
            % error
            elseif violin_index == num_individual_methods + 1
                % set color
                violins(violin_index).FaceColor = [0, 0, 0];
                % extract rms
                y = err_rms(component_index);
                % plot dashed line for rms error
                plot([min(violins(1).XData), max(violins(numel(violins)).XData)], y*[1, 1], '--', 'color', [0, 0, 0, 0.2])
            
            % combined uncertainty methods
            elseif violin_index > num_individual_methods + 1
                % set color
                violins(violin_index).FaceColor = colors(2, :);
                % extract rms value
                y = sigma_rms_comb(violin_index - (num_individual_methods + 1), component_index);
            end

            % adjust violin properties
            color_all{violin_index} = violins(violin_index).FaceColor;
            violins(violin_index).FaceAlpha = 0.25;

            % plot rms value
            plot([min(x), max(x)], y*[1, 1], 'color', [violins(violin_index).FaceColor, violins(violin_index).FaceAlpha+0.2], 'linewidth', 3)

            % add rms as text
            x_text = mean(x) + 0.1 * (max(x) - min(x));        
            y_text = y + 0.25 * max_error_threshold(component_index);
            text(x_text, y_text, num2str(y, '%.2f'), 'fontweight', 'bold', 'color', violins(violin_index).FaceColor); 
        end

        % annotate y axis
        if component_index == 2
            yl = ylabel('Error and Uncertainty (pix.)', 'fontsize', 14);
            set(yl, 'units', 'normalized', 'position', [-0.1 0.45 -1]);
        end

        % adjust limits
        ylim([0 max_error_threshold(component_index)])
        pause(0.1);
        drawnow();

        ax = gca;
        if component_index < num_components
            ax.XAxis.Visible = 'off';
        else
            % turn off x axis line
            ax.XAxis.Axle.Visible = 'off';
            ax.XAxis.TickLength = [0 0];
        end
        % turn off y axis line
        ax.YAxis.Axle.Visible = 'off';
        % ax.YAxis.TickLength = [0, 0];
        ax.YAxis.TickLength = [0.005 0.005];
        title(component_names_plot{component_index});
        % text(-0.15, max_error_threshold(component_index)/2, component_names_plot{component_index}, 'fontsize', 14, 'fontweight', 'bold')

        drawnow();
    end

    % adjust figure position
    set(gcf, 'resize', 'off');
    % set(gcf, 'units', 'inches', 'Position', [244  68  1175  873]/user_screen_resolution)
    % set(gcf, 'units', 'inches', 'Position', [244  400  750  539]/user_screen_resolution)
    set(gcf, 'units', 'inches', 'Position', [244  400  650  435]/user_screen_resolution)
    set(gcf, 'resize', 'off');
    drawnow();
    % set(gca, 'units', 'pix', 'fontsize', 11);     
end

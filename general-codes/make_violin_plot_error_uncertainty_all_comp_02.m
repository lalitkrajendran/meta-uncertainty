function violins = make_violin_plot_error_uncertainty_all_comp_02(err, sigma_indv, sigma_comb, err_rms, sigma_rms_indv, sigma_rms_comb, bins, ...
                                                                    individual_method_array, combination_method_array, ...
                                                                    component_names, colors, user_screen_resolution, max_error_threshold)

    num_components = numel(err);
    num_individual_methods = numel(individual_method_array);
    num_combined_methods = numel(combination_method_array);
    num_methods_total = num_individual_methods + num_combined_methods + 1;
    
    figure
    drawnow();    
    % loop through components
    for comp_index = 1:num_components
        subplot(num_components, 1, comp_index)

        % --------------------------
        % obtain max number of elements in any of the arrays
        % --------------------------
        N_all = [numel(err{comp_index})];
        for method_index = 1:num_individual_methods
            N_all = [N_all, numel(sigma_indv{method_index, comp_index})];
        end

        for method_index = 1:num_combined_methods
            N_all = [N_all, numel(sigma_comb{method_index, comp_index})];
        end
        N_max = max(N_all);

        % --------------------------
        % pad arrays with nans to ensure they are the same size
        % --------------------------
        err{comp_index} = padarray(err{comp_index}, [0, N_max - N_all(1)], NaN, 'post');
        for method_index = 1:num_individual_methods
            sigma_indv{method_index, comp_index} = padarray(sigma_indv{method_index, comp_index}, [0, N_max - N_all(method_index+1)], NaN, 'post');
        end

        for method_index = 1:num_combined_methods
            sigma_comb{method_index, comp_index} = padarray(sigma_comb{method_index, comp_index}, [0, N_max - N_all(method_index+4)], NaN, 'post');
        end

        % make violin plot arrays and colors
        plot_array = nans(N_max, num_individual_methods + num_combined_methods + 1);
        colors_array = nans(num_methods_total, 3);
        for method_index = 1:num_individual_methods
            plot_array(:, method_index) = sigma_indv{method_index, comp_index}';
            colors_array(method_index, :) = colors(1, :);
        end

        plot_array(:, num_individual_methods+1) = abs(err{comp_index})';
        colors_array(num_individual_methods+1, :) = [0, 0, 0];

        for method_index = 1:num_combined_methods
            plot_array(:, method_index + num_individual_methods + 1) = sigma_comb{method_index, comp_index}';
            colors_array(method_index + num_individual_methods + 1, :) = colors(2, :);
        end

        violins = violin_plot_pdf(plot_array, bins{comp_index}, {individual_method_array{:}, 'Error', combination_method_array{:}}, colors_array);

        % face colors
        color_all = cell(1, numel(violins));

        % add lines corresponding to rms
        for violin_index = 1:numel(violins)
            % x co-ordinates of current violin plot
            x = violins(violin_index).XData;

            % individual uncertainty methods
            if violin_index <= num_individual_methods
                violins(violin_index).FaceColor = colors(1, :);
                % violins(violin_index).FaceColor = colors_blue{violin_index-1};
                y = sigma_rms_indv(violin_index, comp_index);
            % error
            elseif violin_index == num_individual_methods + 1
                violins(violin_index).FaceColor = [0, 0, 0];
                y = err_rms(comp_index);
                % x = [violins(violin_index).XData; violins(numel(violins)).XData];
                plot([min(violins(1).XData), max(violins(numel(violins)).XData)], y*[1, 1], '--', 'color', [0, 0, 0, 0.2])
            % combined uncertainty methods
            elseif violin_index > num_individual_methods + 1
                violins(violin_index).FaceColor = colors(2, :);
                % violins(violin_index).FaceColor = colors_red{violin_index-4};
                y = sigma_rms_comb(violin_index - (num_individual_methods + 1), comp_index);
            end

            %% adjust violin properties
            color_all{violin_index} = violins(violin_index).FaceColor;
            violins(violin_index).FaceAlpha = 0.25;

            %% plot rms value
            plot([min(x), max(x)], y*[1, 1], 'color', [violins(violin_index).FaceColor, violins(violin_index).FaceAlpha+0.2], 'linewidth', 3)
        end

        % annotate y axis
        if comp_index == 2
            yl = ylabel('Error/Uncertainty (pix.)', 'fontsize', 16);
            set(yl, 'units', 'normalized', 'position', [-0.05 0.1 -1]);
        end
        
        % adjust limits
        ylim([0 max_error_threshold(comp_index)])
        pause(0.1);
        drawnow();

        ax = gca;
        if comp_index < num_components
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
        % title(component_names{comp_index});
        text(-0.15, max_error_threshold(comp_index)/2, component_names{comp_index}, 'fontsize', 14, 'fontweight', 'bold')
    end
    
    % adjust figure position
    set(gcf, 'resize', 'off');
    set(gcf, 'units', 'inches', 'Position', [244  68  1175  873]/user_screen_resolution)
    set(gcf, 'resize', 'off');
    drawnow();
    % set(gca, 'units', 'pix', 'fontsize', 11);    
end

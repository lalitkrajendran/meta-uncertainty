function violins = make_violin_plot_error_uncertainty_all_comp(err_abs, sigma_indv, sigma_comb, err_rms, sigma_rms_indv, sigma_rms_comb, ...
                                                                individual_method_array, combination_method_array, ...
                                                                component_names, colors, user_screen_resolution, max_error_threshold)
    
    num_components = numel(err_abs);

    figure
    set(gcf, 'resize', 'off');
    drawnow();
    set(gcf, 'units', 'inches', 'Position', [244          68        1175         873]/user_screen_resolution)
    drawnow();

    % loop through components
    for comp_index = 1:num_components
        subplot(num_components, 1, comp_index)
        % obtain max number of elements in any of the arrays
        N_all = [numel(err_abs{comp_index}), numel(sigma_indv{1, comp_index}), numel(sigma_indv{2, comp_index}), numel(sigma_indv{3, comp_index}), ...
                    numel(sigma_comb{1, comp_index}), numel(sigma_comb{2, comp_index}), numel(sigma_comb{3, comp_index})];
        N_max = max(N_all);

        % pad arrays with nans to ensure they are the same size
        err_abs{comp_index} = padarray(err_abs{comp_index}, [0, N_max - N_all(1)], NaN, 'post');
        for method_index = 1:3
            sigma_indv{method_index, comp_index} = padarray(sigma_indv{method_index, comp_index}, [0, N_max - N_all(method_index+1)], NaN, 'post');
        end
        for method_index = 1:3
            sigma_comb{method_index, comp_index} = padarray(sigma_comb{method_index, comp_index}, [0, N_max - N_all(method_index+4)], NaN, 'post');
        end

        % make violin plot
        violins = violinplot([sigma_indv{1, comp_index}', sigma_indv{2, comp_index}', sigma_indv{3, comp_index}', abs(err_abs{comp_index})', ...
                                sigma_comb{1, comp_index}', sigma_comb{2, comp_index}', sigma_comb{3, comp_index}'], ...
                                {individual_method_array{:}, 'Error', combination_method_array{:}} , ...
                                'showdata', false, 'shownotches', false, 'edgecolor', [1, 1, 1], 'ViolinAlpha', 0.5);

        % face colors
        color_all = cell(1, numel(violins));

        % add lines corresponding to rms
        for violin_index = 1:numel(violins)
            % x co-ordinates of current violin plot
            x = violins(violin_index).ViolinPlot.XData;

            % individual uncertainty methods
            if violin_index <= 3
                violins(violin_index).ViolinColor = colors(1, :);
                % violins(violin_index).ViolinColor = colors_blue{violin_index-1};
                y = sigma_rms_indv(violin_index, comp_index);
            % error
            elseif violin_index == 4
                violins(violin_index).ViolinColor = [0, 0, 0];
                y = err_rms(comp_index);
                % x = [violins(violin_index).ViolinPlot.XData; violins(numel(violins)).ViolinPlot.XData];
                plot([min(violins(1).ViolinPlot.XData), max(violins(numel(violins)).ViolinPlot.XData)], y*[1, 1], '--', 'color', [0, 0, 0, 0.2])
            % combined uncertainty methods
            elseif violin_index > 4
                violins(violin_index).ViolinColor = colors(2, :);
                % violins(violin_index).ViolinColor = colors_red{violin_index-4};
                y = sigma_rms_comb(violin_index - 4, comp_index);
            end

            %% adjust violin properties
            color_all{violin_index} = violins(violin_index).ViolinColor;
            violins(violin_index).ViolinAlpha = 0.25;
            violins(violin_index).BoxColor = violins(violin_index).ViolinColor;
            violins(violin_index).BoxWidth = 0.005;

            %% plot rms value
            plot([min(x), max(x)], y*[1, 1], 'color', [violins(violin_index).ViolinColor, violins(violin_index).ViolinAlpha+0.2], 'linewidth', 3)

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

        if comp_index == 1
            text(1.5, max_error_threshold(comp_index)+0.1, 'Individual',  'fontsize', 14, 'fontweight', 'bold', 'color', colors(1, :));
            text(5.5, max_error_threshold(comp_index)+0.1, 'Combined',  'fontsize', 14, 'fontweight', 'bold', 'color', colors(2, :));
        end

    end
    % adjust figure position
    set(gcf, 'resize', 'off');
    drawnow();
    % set(gcf, 'units', 'inches', 'Position', [352   526   895   392]/user_screen_resolution)
    set(gcf, 'units', 'inches', 'Position', [244          68        1175         873]/user_screen_resolution)
    drawnow();
    % set(gca, 'units', 'pix', 'fontsize', 11);    
end

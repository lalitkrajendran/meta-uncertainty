function violins = make_violin_plot_error_uncertainty(err_abs, sigma_indv, sigma_comb, err_rms, sigma_rms_indv, sigma_rms_comb, ...
                                                        individual_method_array, combination_method_array, ...
                                                        colors, user_screen_resolution, max_error_threshold)
        
    % number of methods
    num_individual_methods = numel(individual_method_array);
    num_combined_methods = numel(combination_method_array);

    % obtain max number of elements in any of the arrays
    N_all = nans(1, num_individual_methods + num_combined_methods + 1);    
    N_all(1) = numel(err_abs);
    for i = 1:num_individual_methods
        N_all(1 + i) = numel(sigma_indv{i});
    end
    for i = 1:num_combined_methods
        N_all(1 + num_individual_methods + i) = numel(sigma_comb{i});
    end
    N_max = max(N_all);

    % pad arrays with nans to ensure they are the same size
    err_abs = padarray(err_abs, [0, N_max - N_all(1)], NaN, 'post');
    for method_index = 1:num_individual_methods
        sigma_indv{method_index} = padarray(sigma_indv{method_index}, [0, N_max - N_all(method_index+1)], NaN, 'post');
    end
    for method_index = 1:num_combined_methods
        sigma_comb{method_index} = padarray(sigma_comb{method_index}, [0, N_max - N_all(method_index+1+num_individual_methods)], NaN, 'post');
    end

    % make violin plot
    plot_array = nans(N_max, num_individual_methods + num_combined_methods + 1);
    for i = 1:num_individual_methods
        plot_array(:, i) = sigma_indv{i}';
    end
    plot_array(:, num_individual_methods+1) = err_abs';
    for i = 1:num_combined_methods
        plot_array(:, i + num_individual_methods + 1) = sigma_comb{i}';
    end

    figure
    % set(gcf, 'visible', 'off');
    % violins = violinplot([err_abs', sigma_indv{1}', sigma_indv{2}', sigma_indv{3}', ...
    % violins = violinplot([sigma_indv{1}', sigma_indv{2}', sigma_indv{3}', err_abs', ...
    %                         sigma_comb{1}', sigma_comb{2}', sigma_comb{3}'], ...
    %                         {individual_method_array{:}, 'Error', combination_method_array{:}} , ...
    %                     'showdata', false, 'shownotches', false, 'edgecolor', [1, 1, 1], 'ViolinAlpha', 0.5);

    violins = violinplot(plot_array, {individual_method_array{:}, 'Error', combination_method_array{:}} , ...
                        'showdata', false, 'shownotches', false, 'showmean', true, 'edgecolor', [1, 1, 1], 'ViolinAlpha', 0.5);

    % face colors
    color_all = cell(1, numel(violins));

    % add lines corresponding to rms
    for violin_index = 1:numel(violins)
        % x co-ordinates of current violin plot
        x = violins(violin_index).ViolinPlot.XData;

        % individual uncertainty methods
        if violin_index <= num_individual_methods
            violins(violin_index).ViolinColor = colors(1, :);
            % violins(violin_index).ViolinColor = colors_blue{violin_index-1};
            y = sigma_rms_indv(violin_index);
        % error
        elseif violin_index == num_individual_methods+1
            violins(violin_index).ViolinColor = [0, 0, 0];
            y = err_rms;
            % x = [violins(violin_index).ViolinPlot.XData; violins(numel(violins)).ViolinPlot.XData];
            plot([min(violins(1).ViolinPlot.XData), max(violins(numel(violins)).ViolinPlot.XData)], y*[1, 1], '--', 'color', [0, 0, 0, 0.2])
        % combined uncertainty methods
        elseif violin_index > num_individual_methods + 1
            violins(violin_index).ViolinColor = colors(2, :);
            % violins(violin_index).ViolinColor = colors_red{violin_index-4};
            y = sigma_rms_comb(violin_index - (num_individual_methods + 1));
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
    ylabel('Error/Uncertainty (pix.)', 'fontsize', 16)
    % adjust limits
    ylim([0 max_error_threshold])
    % ylim([0 0.15])

    pause(0.1);
    % turn off x axis line
    ax = gca;
    ax.XAxis.Axle.Visible = 'off';
    ax.XAxis.TickLength = [0 0];
    % turn off y axis line
    ax.YAxis.Axle.Visible = 'off';
    % ax.YAxis.TickLength = [0, 0];
    ax.YAxis.TickLength = [0.005 0.005];

    % adjust figure position
    set(gcf, 'resize', 'off');
    drawnow();
    set(gcf, 'units', 'inches', 'Position', [352   526   895   392]/user_screen_resolution)
    drawnow();
    % set(gca, 'units', 'pix', 'fontsize', 11);    
end

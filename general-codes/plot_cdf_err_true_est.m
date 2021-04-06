function plot_cdf_err_true_est(bins, cdf_err, cdf_err_est_individual, cdf_err_est_combined, max_error_threshold, uncertainty_method_array, combination_method_names, ...
                            colors, line_symbols, user_screen_resolution)
    num_uncertainty_methods = numel(uncertainty_method_array);
    num_combination_methods = numel(combination_method_names);
                        
    figure
    ax1 = subplot(2, 1, 1);
    for uncertainty_method_index = 1:num_uncertainty_methods    
        plot(bins(1:end-1), cdf_err_est_individual(uncertainty_method_index, :), line_symbols{uncertainty_method_index}, 'color', colors(1, :))
        hold on
    end
    plot(bins(1:end-1), cdf_err, 'k')
    xlim([0 max_error_threshold])
    ylim([0 1])
    box off
    set(gca, 'xticklabel', '');
    legend({uncertainty_method_array{:}, 'Error'}, 'location', 'southeast')

    ax2 = subplot(2, 1, 2);
    for combination_method_index = 1:num_combination_methods  
        plot(bins(1:end-1), cdf_err_est_combined(combination_method_index, :), line_symbols{combination_method_index}, 'color', colors(2, :))
        hold on
    end
    plot(bins(1:end-1), cdf_err, 'k')
    xlim([0 max_error_threshold])
    ylim([0 1])
    box off
    xlabel('(pix.)', 'fontsize', 16)
    legend({combination_method_names{:}, 'Error'}, 'location', 'southeast')

    annotation('textarrow',[0.5 0.5],[0.5 0.5],'string','Cumulative Density Function (CDF)', ...
    'HeadStyle','none','LineStyle', 'none', 'TextRotation',90,'Position',[.03 .85 0 0],'Fontsize',16); %,'FontWeight','bold');
    set(gcf, 'resize', 'off');
    set(gcf, 'units', 'inches', 'position', [680   402   500   600]/user_screen_resolution);
    set(gca, 'units', 'pix', 'fontsize', 11);

end
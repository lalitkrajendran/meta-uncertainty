function plot_pdf_uncertainties_individual_resampled(pdf_unc_individual, pdf_unc_resampled, rms_unc_individual, rms_unc_resampled, bins, method_names, user_screen_resolution)
    % number of methods
    num_methods = numel(method_names);
    % get line colors
    colors = lines(1);
    % calculate y limit for the plots
    y_max = max([pdf_unc_individual(:); pdf_unc_resampled(:)]) * 1.1;
    
    % make plots
    figure
    for method_index = 1:num_methods
        subplot(num_methods, 1, method_index)
        l1 = plot(bins(1:end-1), pdf_unc_individual(method_index, :), 'color', colors(1, :));
        hold on
        plot(rms_unc_individual(method_index) * [1, 1], [0, y_max], 'color', colors(1, :))
        l2 = plot(bins(1:end-1), pdf_unc_resampled(method_index, :), '--', 'color', colors(1, :));
        plot(rms_unc_resampled(method_index) * [1, 1], [0, y_max], '--', 'color', colors(1, :))

        xlim([0 bins(end)])
        ylim([0 y_max])
        box off
        if method_index == num_methods
            xlabel('(pix.)')
        else
            set(gca, 'xticklabel', []);
        end
        title(method_names{method_index})
        if method_index == 1
            legend([l1, l2], 'Individual', 'Resampled')
        end
    end

    % adjsut plot
    set(gcf, 'resize', 'off');
    set(gcf, 'units', 'inches', 'position', [680   370   530   480]/user_screen_resolution);
    set(gca, 'units', 'pix', 'fontsize', 11);

end
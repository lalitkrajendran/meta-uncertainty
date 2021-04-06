function lgd = plot_binwise_rms(err_rms_binned, unc_indiv_rms_binned, unc_comb_rms_binned, num_indices, bins, ...
                                individual_method_names, resampling_method_names, symbols, colors, user_screen_resolution)
    
    % calculate number of individual methods
    num_individual_methods = size(unc_indiv_rms_binned, 2);
    % calculate number of combined methods
    num_combined_methods = size(unc_comb_rms_binned, 2);
    % calculate number of bins
    num_bins = numel(bins);

    figure
    subplot(2, 1, 1)
    plot(bins, err_rms_binned, 'k')
    hold on
    for method_index = 1:num_individual_methods
        plot(bins, unc_indiv_rms_binned(:, method_index), symbols{method_index}, 'color', colors(1, :))
    end
    plot(bins, unc_comb_rms_binned, symbols{1}, 'color', colors(2, :))
    box off
    xlim([-0.5 0.5])
    ylim([0 0.2])
    set(gca, 'xticklabel', [])
    % xlabel('Fractional Displacement (pix.)')
    ylabel({'RMS Error and'; 'Uncertainty (pix.)'}, 'fontsize', 10)
    lgd = legend({'Error', individual_method_names{:}, resampling_method_names{:}}, 'location', 'northoutside', 'orientation', 'horizontal');
    lgd.FontSize = 9;

    subplot(2, 1, 2)
    bar(bins, num_indices/sum(num_indices) * 100)
    box off
    xlim([bins(1) bins(end)])
    ylim([0 100])
    % ylim([0 num_valid_trials_consolidated])
    xlabel('Fractional Displacement (pix.)')
    ylabel('Bin Count (%)', 'fontsize', 10)

    set(gcf, 'resize', 'off');
    set(gcf, 'units', 'inches', 'position', [274   304   450   477]/user_screen_resolution);
    set(gcf, 'resize', 'off');
    drawnow();

    lgd.Position(1) = 0;
    lgd.Position(2) = 0.95;

end
    
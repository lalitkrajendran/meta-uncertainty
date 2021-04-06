function plot_snr_ratio_violin(snr_resampling_trials, snr_sub_trials, ppr, num_resampling_trials, snr_method_names, user_screen_resolution)

    num_snr_methods = numel(snr_method_names);        
    num_ppr = numel(ppr);

    snr_ratio = nans(num_snr_methods, num_ppr, num_resampling_trials);    
    ppr_name_array = cell(1, num_ppr);

    % ==========================
    % loop through uncertainty methods
    % ==========================    
    for method_index = 1:num_snr_methods
        % method name
        method_name = snr_method_names{method_index};
        % extract original uncertainties
        snr_sub = snr_sub_trials.(method_name);
        % ================================================
        %% loop through particle removal percentages
        % ================================================
        for ppr_index = 1:num_ppr
            snr_ratio(method_index, ppr_index, :) = snr_resampling_trials{ppr_index}.(method_name) / snr_sub;
            % name for this case
            if method_index == 1
                ppr_name_array{ppr_index} = num2str(round(ppr(ppr_index)), '%d');
            end
        end
    end

    % bins = linspace(0.5, 1.5, 10);
    bins = linspace(0, 2, 10);
    figure
    for method_index = 1:num_snr_methods
        % method name
        method_name = (snr_method_names{method_index});
        subplot(1, num_snr_methods, method_index)
        
        % violinplot_single_symm(snr_ratio, particle_remove_index, 0.3, lines(1));
        % violins = violin_plot_pdf(squeeze(snr_ratio(method_index, :, :)), ppr_name_array, 'violincolor', colors(1, :), 'showdata', false, 'shownotches', true, 'showmean', true, ...
        %                     'showrms', true, 'edgecolor', [1, 1, 1], 'ViolinAlpha', 0.5);
        plot_array = squeeze(snr_ratio(method_index, :, :))';
        violins = violin_plot_pdf(plot_array, bins, ppr_name_array, repmat(lines(1), num_ppr, 1));

        ylim([0, 2])
        
        if method_index == 1
            ylabel('SNR Ratio')
        else            
            set(gca, 'yticklabels', [])
        end
        
        xlabel('Particle Addition %')
        title(method_name)
        % % --------------------------
        % % annotate figure
        % % --------------------------                    
        % xlim([1 num_ppr+0.5])
        % ylim([0 2])

        % if method_index == num_snr_methods
        %     if resampling_method_index <= 2
        %         xlabel('Particle Removal %');
        %     else
        %         xlabel('Particle Addition %');
        %     end
        % else
        %     set(gca, 'xticklabel', []);
        % end

        % if resampling_method_index == 1
        %     ylabel('Uncertainty Ratio')                    
        % else
        %     set(gca, 'yticklabel', []);
        % end

        % if method_index == 1
        %     title(resampling_method_names{resampling_method_index})
        % end                    

        % if resampling_method_index == 1
        %     annotation('textbox', [0.02, 0.85 - 0.3 * (method_index- 1), 0, 0], 'string', upper(method_name), 'fontsize', 18, 'fontweight', 'bold')
        % end

    end
    set(gcf, 'resize', 'off');
    set(gcf, 'units', 'inches', 'position', [233   442   520   278]/user_screen_resolution)
    set(gcf, 'resize', 'off');
    drawnow();

end
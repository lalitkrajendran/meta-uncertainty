function plot_unc_ratio_violin(unc_resampling_trials, unc_sub_trials, ppr, num_resampling_trials, individual_method_names, user_screen_resolution)

    num_individual_methods = numel(individual_method_names);        
    num_ppr = numel(ppr);

    unc_ratio = nans(num_individual_methods, num_ppr, num_resampling_trials);
    component_name = 'x';
    ppr_name_array = cell(1, num_ppr);
    
    colors = lines(3);

    % ==========================
    % loop through uncertainty methods
    % ==========================    
    for individual_method_index = 1:num_individual_methods
        % method name
        method_name = lower(individual_method_names{individual_method_index});
        % extract original uncertainties
        unc_sub = unc_sub_trials{individual_method_index}.(component_name);
        % ================================================
        %% loop through particle removal percentages
        % ================================================
        for ppr_index = 1:num_ppr
            unc_ratio(individual_method_index, ppr_index, :) = unc_resampling_trials{ppr_index}.([method_name component_name]) / unc_sub;
            % name for this case
            if individual_method_index == 1
                ppr_name_array{ppr_index} = num2str(round(ppr(ppr_index)), '%d');
            end
        end
    end

    % bins = linspace(0.5, 1.5, 10);
    bins = linspace(0, 2, 10);
    figure
    for individual_method_index = 1:num_individual_methods
        % method name
        method_name = (individual_method_names{individual_method_index});
        subplot(1, num_individual_methods, individual_method_index)
        
        % violinplot_single_symm(unc_ratio, particle_remove_index, 0.3, lines(1));
        % violins = violin_plot_pdf(squeeze(unc_ratio(individual_method_index, :, :)), ppr_name_array, 'violincolor', colors(1, :), 'showdata', false, 'shownotches', true, 'showmean', true, ...
        %                     'showrms', true, 'edgecolor', [1, 1, 1], 'ViolinAlpha', 0.5);
        plot_array = squeeze(unc_ratio(individual_method_index, :, :))';
        violins = violin_plot_pdf(plot_array, bins, ppr_name_array, repmat(lines(1), num_ppr, 1));

        ylim([0, 2])
        
        if individual_method_index == 1
            ylabel('Uncertainty Ratio')
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

        % if individual_method_index == num_individual_methods
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

        % if individual_method_index == 1
        %     title(resampling_method_names{resampling_method_index})
        % end                    

        % if resampling_method_index == 1
        %     annotation('textbox', [0.02, 0.85 - 0.3 * (individual_method_index- 1), 0, 0], 'string', upper(method_name), 'fontsize', 18, 'fontweight', 'bold')
        % end

    end


    % plot invisible markers for legend
    l1 = plot(NaN, NaN, '^', 'color', colors(1, :));
    l2 = plot(NaN, NaN, 'o', 'color', colors(1, :));
    l3 = plot(NaN, NaN, 'v', 'color', colors(1, :));
    lgd = legend([l1, l2, l3], 'Q1', 'Median', 'Q3', 'location', 'northoutside', 'Orientation', 'horizontal');

    set(gcf, 'resize', 'off');
    % set(gcf, 'units', 'inches', 'position', [233   442   781   278]/user_screen_resolution)
    set(gcf, 'units', 'inches', 'position', [233   442   781   300]/user_screen_resolution)
    set(gcf, 'resize', 'off');
    drawnow();

    % set legend position
    lgd.Position(1:2) = [0.3, 0.9];
    set_subplots_position_y(gcf, 0.15, 0.65);


    drawnow();

end
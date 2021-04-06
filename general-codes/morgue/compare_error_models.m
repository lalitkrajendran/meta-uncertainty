function compare_error_models(error_models, err, err_est_indiv, err_est_comb, individual_method_array, resampling_method_names_plot_short, bins, user_screen_resolution)

        num_models = numel(error_models);
        num_individual_methods = numel(individual_method_array);
        num_resampling_methods = numel(resampling_method_names_plot_short);

        colors = lines(3);

        p1 = histcounts(abs(err), bins, 'normalization', 'pdf');
        h = figure;
        method_type_names = {'Indiv'; 'Comb'}; 
        for method_type = [1, 2]            
            if method_type == 1
                num_methods = num_individual_methods;
                err_est = err_est_indiv;
                method_names = individual_method_array;
            else
                num_methods = num_resampling_methods;
                err_est = err_est_comb;
                method_names = resampling_method_names_plot_short;
            end

            for method_index = 1:num_methods
                method_name = method_names{method_index};
                subplot_index = (method_type - 1) * num_methods + method_index;
                subplot(2, num_methods, subplot_index)
                plot(bins(1:end-1), p1, 'color', 'k')
                hold on
                
                for model_index = 1:num_models
                    model = error_models{model_index};
                    p = histcounts(abs(err_est{method_index, model_index}), bins, 'normalization', 'pdf');
                    plot(bins(1:end-1), p, 'color', colors(model_index, :))
                end

                box off
                if method_index > 1
                    set(gca, 'yticklabels', [])
                end
                if method_type == 1
                    set(gca, 'xticklabels', [])
                else
                    xlabel('Error Magnitude (pix.)')
                end
                title(method_name)

                if method_index == 1
                    annotation('textbox', [0.01, 0.8 - 0.5 * (method_type - 1), 0, 0], 'string', method_type_names{method_type}, 'fontsize', 14, 'fontweight', 'bold')    
                end
            end
            
        end

        figure(h);
        set_common_limits_subplot(gcf, 'x', 'manual', [0, bins(end)]);
        set_common_limits_subplot(gcf, 'y', 'auto');

        subplot(2, num_methods, 1)
        legend('Error', error_models{:})

        set(gcf, 'resize', 'off');
        set(gcf, 'units', 'inches', 'position', [144   448   963   550]/user_screen_resolution);
        set(gcf, 'resize', 'off');
        drawnow();

end
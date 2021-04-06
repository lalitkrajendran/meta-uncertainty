function plot_consolidated_qq_v3(results, dataset_name_array, individual_method_array, combined_method_array, combined_method_index_plot, max_error_threshold, symbols, user_screen_resolution)
    
    % extract number of window resolutions
    num_window_resolution = size(results, 1);
    % extract number of datasets
    num_datasets = size(results, 2);
    % extract number of individual methods
    num_individual_methods = numel(individual_method_array);
    % extract number of combined methods
    num_combined_methods = numel(combined_method_array);
    % calculate total number of methods
    num_total_methods = num_individual_methods + num_combined_methods + 1;
    % generate color array
    colors = lines(3);
    % plot marker size
    marker_size = 4;
    % number of data points to skip for the qq plot
    num_skip = 24;

    % loop through window resolutions
    for window_resolution_index = 1:num_window_resolution
        % loop through datasets
        for dataset_index = 1:num_datasets
            % name of the current dataset
            dataset_name = dataset_name_array{dataset_index};
            
            % =============================
            %% qq plot of error and uncertainty histograms
            % =============================
            subplot_index = (window_resolution_index - 1) * num_datasets + dataset_index;
            subplot(2, num_datasets, subplot_index);
            ax{subplot_index} = gca;
            plot([-1 1] * max_error_threshold, [-1 1] * max_error_threshold, 'k');
            hold on

            x = results{window_resolution_index, dataset_index}.err_all;
            x = x(isfinite(x));
            legend_string = cell(1, num_individual_methods + num_combined_methods);
            % individual methods
            qq_indiv = [];
            for method_index = 1:num_individual_methods
                y = results{window_resolution_index, dataset_index}.err_est_indiv(:, method_index);
                l = make_qq_plot(x, y, num_skip, colors(1, :), symbols{method_index}, marker_size);
                
                qq_indiv = [qq_indiv, l(1)];
            end
            % combined methods
            qq_comb = [];
            qq_rms_sigma_comb = nans(1, num_combined_methods);
            for method_index = combined_method_index_plot
                y = results{window_resolution_index, dataset_index}.err_est_comb(:, method_index);
                % l = make_qq_plot(x, y, num_skip, colors(2, :), symbols{method_index}, marker_size);
                l = make_qq_plot(x, y, num_skip, colors(2, :), symbols{1}, marker_size);
                
                qq_comb = [qq_comb, l(1)];
            end

            box off
            axis equal

            axis([-max_error_threshold, max_error_threshold, -max_error_threshold, max_error_threshold])
            set(gca, 'YAxisLocation', 'origin')
            set(gca, 'XAxisLocation', 'origin')
        
            % add dataset name
            if window_resolution_index == 1
                title(dataset_name);
                set(gca, 'Position', [ax{subplot_index}.Position(1)-0.03, 0.5, 0.15, 0.35]);
                if dataset_index == 3
                    % add legend
                    lgd = legend([qq_indiv, qq_comb], ...
                    {individual_method_array{:}, combined_method_array{:}}, 'location', 'northoutside', 'orientation', 'horizontal');
                    lgd.Position = [0.25    0.92    0.5490    0.0416];
                    lgd.FontSize = 10;
                end
            elseif window_resolution_index == 2
                set(gca, 'Position', [ax{subplot_index}.Position(1)-0.03, 0.1, 0.15, 0.35]);
            end
            
            % if window_resolution_index == 2 && dataset_index == 1
            if dataset_index == 1
                % xlabel('True Error (pix.)', 'fontsize', 10);
                % ylabel('Estimated Error (pix.)', 'fontsize', 10);
                % xl = xlabel('True Error (pix.)', 'fontsize', 10); %16)
                xl = xlabel('\epsilon_{True} (pix.)', 'fontsize', 10); %16)
                set(xl, 'HorizontalAlignment', 'center');
                % yl = ylabel('Estimated Error (pix.)', 'fontsize', 10); %16)
                yl = ylabel('\epsilon_{Est} (pix.)', 'fontsize', 10); %16)
                set(yl, 'HorizontalAlignment', 'right');            

                text(-0.4, 0, ['WS ' num2str(window_resolution_index)], 'fontsize', 14, 'fontweight', 'bold')
            else
                xlabel('');
                ylabel('');
                ax{subplot_index}.XAxis.TickLabels = '';
                ax{subplot_index}.YAxis.TickLabels = '';
                pause(0.1);
            end        
        end
    end

    set(gcf, 'resize', 'off');
    pause(0.1);
    set(gcf, 'units', 'inches', 'position', [185   257   1300   650]/user_screen_resolution);
    set(gcf, 'resize', 'off');
    drawnow();

end

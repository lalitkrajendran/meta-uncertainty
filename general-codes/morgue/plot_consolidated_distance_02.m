function plot_consolidated_distance_02(results, dataset_name_array, individual_method_array, combined_method_array)
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

    % categories
    histogram_distance_categories = categorical({individual_method_array{:}, combined_method_array{:}});
    % sort in the desired order
    histogram_distance_categories = reordercats(histogram_distance_categories, ...
                                    {individual_method_array{:}, combined_method_array{:}});

    % loop through window resolutions
    for window_resolution_index = 1:num_window_resolution
        % loop through datasets
        for dataset_index = 1:num_datasets
            % name of the current dataset
            dataset_name = dataset_name_array{dataset_index};

            % ===================================
            % plot total variation distance
            % ===================================
            distance_method_index = 1;

            subplot_index = (window_resolution_index - 1) * num_datasets + dataset_index;
            subplot(2, num_datasets, subplot_index);
            ax{subplot_index} = gca;
            % create array of distances
            Y = [results{window_resolution_index, dataset_index}.d_err_est_indiv, results{window_resolution_index, dataset_index}.d_err_est_comb];
            X = histogram_distance_categories;

            % find max value
            [maxval, maxloc] = max(Y);
            % find min value
            [minval, minloc] = min(Y);
            
            % normalize distances to the maximum
            Y = Y/maxval * 100;
            % make bar plot
            b = bar(X,Y); %,'FaceColor','flat');
            hold on
            
            % plot minimum distance level
            plot(X, minval/maxval * 100 * ones(1, numel(X)), '--', 'color', [0, 0, 0, 0.2], 'linewidth', 3.0);

            % adjust bar properties
            b.FaceAlpha = 0.7;
            b.BarWidth = 0.5;
            b.FaceColor = 'flat';
            b.EdgeColor = [1, 1, 1];
            b.ShowBaseLine = 'off';
            % adjust bar color
            for i = 1:numel(X)
                if i <= 3
                    b.CData(i, :) = colors(1, :);
                else
                    b.CData(i, :) = colors(2, :);
                end
            end
            % % add text on top of each bar
            % text(1:length(Y),Y,num2str(Y', '%.1f'),'vert','bottom','horiz','center');
            % set axis limit
            ylim([0 100])
            
            % remove box
            box off
            % annotate axis
            ylabel('(%)')
            
            % remove y axis
            % set(gca, 'ycolor', 'none');
            set(gca,'TickLength',[0 0])
            
            % turn off x axis
            % ax = gca;
            ax{subplot_index}.XAxis.Axle.Visible = 'off';
            % remove ticks
            xtickangle(0)
            
            ax{subplot_index}.Position(1) = ax{subplot_index}.Position(1) - 0.075; 
            ax{subplot_index}.Position(1) = ax{subplot_index}.Position(1) * 1.1; 
            ax{subplot_index}.Position(3) = ax{subplot_index}.Position(3) + 0.05;

            if window_resolution_index == 1
                title(dataset_name_array{dataset_index});
            end
            
            % if ~ (window_resolution_index == 2 && dataset_index == 1)
            if ~ (dataset_index == 1)
                % xlabel('');
                set(gca, 'ycolor', 'none');
                ax{subplot_index}.XAxis.TickLabels = '';
                ax{subplot_index}.YAxis.TickLabels = '';
                pause(0.1);
            end               
        end
    end

    pause(0.1);
    set(gcf, 'resize', 'off');
    set(gcf, 'Position', [29         246        1204         427]);

end


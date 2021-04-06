function plot_consolidated_violins_v5(results, dataset_name_array, individual_method_array, combined_method_array, combined_method_index_plot, max_error_threshold, user_screen_resolution)
% Function to plot violins from all datasets and window resolutions
%
% INPUTS:
% results: results to be plotted (cell array)
% dataset_name_array: names of the datasets
% individual_method_array: individual method names
% 
% OUTPUTS:
% None
%
% AUTHOR:
% Lalit Rajendran (lrajendr@purdue.edu)

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

    % loop through datasets
    for dataset_index = 1:num_datasets
        % name of the current dataset
        dataset_name = dataset_name_array{dataset_index};
        % create axes for the current dataset
        ax = subplot(num_datasets, 1, dataset_index);
        hold on
        % set subplot position
        set(ax, 'Position', [0.15, 0.85 - 0.2 * (dataset_index - 1),  0.65, 0.12]);

        % initialize array to hold violin colors
        violin_colors = cell(1, num_total_methods);
        % initialize array to hold violin plots        
        violins = cell(1, num_total_methods);
        % plot violins and rms values
        for violin_index = 1:num_total_methods
            % error
            if violin_index == 4
                violins{violin_index} = violinplot_single_asymm(abs(results{1, dataset_index}.err_all), ...
                                            abs(results{2, dataset_index}.err_all), violin_index, 0.3, [0, 0, 0]);

                y1 = results{1, dataset_index}.err_rms;
                y2 = results{2, dataset_index}.err_rms;                
            % individual uncertainty methods
            elseif violin_index <= num_individual_methods
                violins{violin_index} = violinplot_single_asymm(results{1, dataset_index}.unc_indiv_all(:, violin_index), ...
                                        results{2, dataset_index}.unc_indiv_all(:, violin_index), violin_index, 0.3, colors(1, :));

                % plot([violin_index, violin_index], [0, max_error_threshold], 'color', [0, 0, 0, 0.5], 'linewidth', 0.5);
                y1 = results{1, dataset_index}.unc_indiv_rms(violin_index);
                y2 = results{2, dataset_index}.unc_indiv_rms(violin_index);
            % combined uncertainty methods
            elseif violin_index > 4
                % violins{violin_index} = violinplot_single_asymm(results{1, dataset_index}.unc_comb_all(:, violin_index-4), ...
                %                         results{2, dataset_index}.unc_comb_all(:, violin_index-4), violin_index, 0.3, colors(2, :));
                violins{violin_index} = violinplot_single_asymm(results{1, dataset_index}.unc_comb_all(:, combined_method_index_plot), ...
                                        results{2, dataset_index}.unc_comb_all(:, combined_method_index_plot), violin_index, 0.3, colors(2, :));

                % plot([violin_index, violin_index], [0, max_error_threshold], 'color', [0, 0, 0, 0.5], 'linewidth', 0.5);
                y1 = results{1, dataset_index}.unc_comb_rms(combined_method_index_plot);
                y2 = results{2, dataset_index}.unc_comb_rms(combined_method_index_plot);
            end

            %% adjust violin properties
            violin_colors{violin_index} = violins{violin_index}(1).FaceColor;
            
            % extract number of plotted points in each half
            num_points = numel(violins{violin_index}(1).XData)/2;
            %% extract plotted points
            x = [violins{violin_index}(1).XData(num_points+1:end); violins{violin_index}(2).XData(1:num_points)];
            y = [violins{violin_index}(1).YData(num_points+1:end); violins{violin_index}(2).YData(1:num_points)];

            %% plot dividing line
            plot([violin_index, violin_index], [min(y), max(y)], 'color', violin_colors{violin_index}, 'linewidth', 0.5);

            %% plot rms value
            % plot([min(x), max(x)], y*[1, 1], 'color', [violins{violin_index}.FaceColor, violins{violin_index}.FaceAlpha+0.2], 'linewidth', 3)
            plot([min(x), violin_index], y1*[1, 1], 'color', [violins{violin_index}(1).FaceColor, violins{violin_index}(1).FaceAlpha+0.2], 'linewidth', 1.5)
            plot([violin_index, max(x)], y2*[1, 1], 'color', [violins{violin_index}(2).FaceColor, violins{violin_index}(2).FaceAlpha+0.2], 'linewidth', 1.5)

            % plot rms error
            plot([min(x), violin_index], results{1, dataset_index}.err_rms*[1, 1], 'color', [0, 0, 0, violins{violin_index}(1).FaceAlpha+0.2], 'linewidth', 2)
            plot([violin_index, max(x)], results{2, dataset_index}.err_rms*[1, 1], 'color', [0, 0, 0, violins{violin_index}(2).FaceAlpha+0.2], 'linewidth', 2)
        end

        % adjust limits
        ylim([0 max_error_threshold])
        
        pause(0.1);

        % turn off x axis line
        ax.XAxis.Axle.Visible = 'off';
        ax.XAxis.TickLength = [0 0];
        if dataset_index < num_datasets
            ax.XAxis.TickLabel = [];
        end

        % turn off y axis line
        ax.YAxis.Axle.Visible = 'off';
        % ax.YAxis.TickLength = [0, 0];
        ax.YAxis.TickLength = [0.005 0.005];

        % add dataset name
        annotation('textbox', [0.85, 0.9 - 0.2 * (dataset_index - 1), 0, 0], 'string', dataset_name_array{dataset_index}, 'fontsize', 14, 'fontweight', 'bold')    
    end

    % xlabel
    xl = set(gca, 'xticklabel', {individual_method_array{:}, 'Error', combined_method_array{:}});

    % annotate y axis
    yl = ylabel('Error and Uncertainty (pix.)', 'fontsize', 16);
    yl.Position = [0.2, 0.65, -1];
    yl.FontWeight = 'bold';
    ax.XAxis.FontSize = 14;
    ax.XAxis.FontWeight = 'bold';

    set(gcf, 'resize', 'off');
    pause(0.1);
    set(gcf, 'units', 'inches', 'position', [332    54   855   690]/user_screen_resolution);
    set(gcf, 'resize', 'off');
    drawnow();
end
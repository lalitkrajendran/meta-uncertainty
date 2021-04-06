function plot_pdf_weights(bins_weights, pdf_w, uncertainty_method_array, combination_method_array, colors, line_symbols, user_screen_resolution)
    % calculate numnber of individual methods
    num_uncertainty_methods = numel(uncertainty_method_array);
    % calculate numnber of combination methods
    num_combination_methods = numel(combination_method_array);
    
    
    figure
    % loop through combination methods
    for combination_method_index = 1:num_combination_methods
        % create subplot for the current method
        % ax(combination_method_index) = subplot(num_combination_methods, 1, combination_method_index);
        % % adjust subplot position
        % set(ax(combination_method_index), 'Position', [0.15, 0.7 - 0.3 * (combination_method_index - 1),  0.65, 0.25]);
        
        % plot weights for unweighted scheme
        if combination_method_index == 1 && strcmpi(combination_method_array{combination_method_index}, 'unwt')
            y_lim = [0, 10];
            plot([1/3, 1/3], y_lim, 'color', colors(1, :))
        % other schemes
        else
            % plot histograms for each uncertainty method
            for uncertainty_method_index = 1:num_uncertainty_methods
                plot(bins_weights(1:end-1), pdf_w{combination_method_index,  uncertainty_method_index}, line_symbols{uncertainty_method_index}, 'color', colors(1, :));
                hold on
            end
            % get y limits
            % y_lim = get(ax(combination_method_index), 'ylim');
            y_lim = get(gca, 'ylim');
            % plot line corresponding to 1/3            
            plot([1/3, 1/3], y_lim, '--', 'color', [0, 0, 0, 0.5])
        end
        % set(gca, 'Position', [0.13, 0.7093, 0.5750, 0.2157])
        box off
        % if combination_method_index < num_combination_methods
        %     set(gca, 'xticklabel', []);
        % else
        %     xl = xlabel('Weights');        
        % end
        xlim([0 1])
        
        % title(combination_method_names{combination_method_index});    
        % annotation('textbox', [0.6, 0.85 - 0.3 * (combination_method_index - 1), 0, 0], 'string', combination_method_names{combination_method_index}, 'fontsize', 16)
        method_name = convert_string_to_sentence_case(combination_method_array{combination_method_index});
        % annotation('textbox', [0.75, 0.85 - 0.3 * (combination_method_index - 1), 0, 0], 'string', method_name, 'fontsize', 14, 'fontweight', 'bold')    
    end

    % axes(ax(combination_method_index))
    lgd = legend(uncertainty_method_array, 'location', 'northoutside', 'orientation', 'horizontal'); %, 'fontsize', 14);

    % lgd.Position = [0.3572 0.9573 0.3199 0.0284];
    % xl.FontSize = 16;
    % % xl.Position = [0.4742  -1   -1.0000];
    % xl.Position = [0.4742  -0.8   -1.0000];

    % annotation('textarrow',[0.5 0.5],[0.5 0.5],'string','Probability Density Function (PDF)', ...
    % 'HeadStyle','none','LineStyle', 'none', 'TextRotation', 90, 'Position', [.04 0.9 0 0], 'FontSize', 14, 'FontWeight', 'bold');
    ylabel('PDF')
    xlabel('Weights');
    set(gcf, 'resize', 'off');
    % set(gcf, 'units', 'inches', 'position', [680   370   530   480]/user_screen_resolution);
    % set(gcf, 'units', 'inches', 'position', [680   370   530   530]/user_screen_resolution);
    set(gcf, 'units', 'inches', 'position', [680  668  515  300]/user_screen_resolution);
    set(gcf, 'resize', 'off');
    drawnow();
    
    set(gca, 'units', 'pix', 'fontsize', 11);

end
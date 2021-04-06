function h = make_histogram_plot(X, Y, colors)
    % normalize distances to the maximum
    % Y = Y/maxval * 100;
    % Y = Y * 100;

    % find max value
    [maxval, maxloc] = max(Y);
    % find min value
    [minval, minloc] = min(Y);
    
    % make bar plot
    b = bar(X,Y); %,'FaceColor','flat');
    hold on
    
    % plot minimum distance level
    % plot(X, minval/maxval * 100 * ones(1, numel(X)), '--', 'color', [0, 0, 0, 0.2], 'linewidth', 3.0);
    plot(X, minval * ones(1, numel(X)), '--', 'color', [0, 0, 0, 0.2], 'linewidth', 3.0);

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
    % add text on top of each bar
    text(1:length(Y),Y,num2str(Y', '%.1f'),'vert','bottom','horiz','center');
    % set axis limit
    % ylim([0 40])
    % ylim([0 100])
    
    % remove box
    box off
    % annotate axis
    ylabel('(%)')
    
    % remove y axis
    set(gca, 'ycolor', 'none');
    set(gca,'TickLength',[0 0])
    
    % turn off x axis
    ax = gca;
    ax.XAxis.Axle.Visible = 'off';
    % remove ticks
    xtickangle(0)
    
    % adjust figure size
    set(gcf, 'resize', 'off');
    pause(0.1);
    set(gcf, 'Position', [100, 300, 600, 375])
    
    % add title
    % t(distance_method_index) = title(convert_string_to_sentence_case(strrep(current_distance_method, '_', ' ')));
    % adjust title position
    % t(distance_method_index).Position(2) = 110;
end
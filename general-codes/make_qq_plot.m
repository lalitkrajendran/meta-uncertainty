function l = make_qq_plot(x, y, num_skip, color, symbol, marker_size)
    % make plot
    % l = qqplot(x(1:num_skip:end), y(1:num_skip:end));
    l = qqplot(x, y, 1:2:100);
    % adjust plot
    % set(l(1), 'marker', 'none'); 
    % set(l(1), 'linestyle', line_symbols{method_index});
    % set(l(1), 'color', colors(1,:));

    set(l(1), 'marker', symbol);
    set(l(1), 'markeredgecolor', color);
    set(l(1), 'linewidth', 0.2);
    % set(l(1), 'markersize', marker_size);
    set(l(1), 'markersize', marker_size);
    
    set(l(2), 'visible', 'off');
    set(l(3), 'visible', 'off');
    
end
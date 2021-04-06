function qq_plot_err_true_est_old(err, err_est_indv, err_est_comb, individual_method_names, combined_method_names, ...
    max_error_threshold, num_skip, colors, symbols, user_screen_resolution)
% Function to make a quantile-quantile plot of the true and estimated error distributions
%
% INPUTS:
% err: true error
% err_est_indv: estimated error from individual uncertainties
% err_est_comb: estimated error from combined methods
% individual_method_names: names of individual uncertainty methods
% combined_method_names: names of combined uncertainty methods
% max_error_threshold: maximum allowable error
% num_skip: number of data points to skip
% colors: marker colors for the individual and combined methods
% symbls: markers
% user_screen_resolution: dpi
%
% OUTPUTS:
% None
%
% AUTHOR:
% Lalit Rajendran (lrajendr@purdue.edu)

    % number of individual methods
    num_individual_methods = numel(individual_method_names);
    % number of combined methods
    num_combined_methods = numel(combined_method_names);

    % plot 1:1 line
    plot([-max_error_threshold max_error_threshold], [-max_error_threshold max_error_threshold], 'k');
    % plot([0 max_error_threshold], [0 max_error_threshold], 'k');
    hold on

    % extract error values
    x = abs(err(:, 1:1000))';
    % only retain finite error values
    x = x(isfinite(x));
    % array to hold plot lines
    l_all = [];
    % plot marker size
    marker_size = 4;

    % individual methods
    qq_sigma = [];
    qq_rms_sigma = nans(1, num_individual_methods);
    % num_skip = 10;
    marker_size = 4;
    legend_string = cell(1, num_individual_methods + num_combined_methods);
    for method_index = 1:num_individual_methods
        y = abs(err_est_indv{method_index}(1:1000, :))';
        % y = err_est_indv{method_index}';

        % make plot
        l = qqplot(x(1:num_skip:end), y(1:num_skip:end));
        % adjust plot
        % set(l(1), 'marker', 'none'); 
        % set(l(1), 'linestyle', line_symbols{method_index});
        % set(l(1), 'color', colors(1,:));

        set(l(1), 'marker', symbols{method_index});
        set(l(1), 'markeredgecolor', colors(1, :));
        set(l(1), 'linewidth', 0.2);
        % set(l(1), 'markersize', marker_size);
        set(l(1), 'markersize', marker_size);

        set(l(2), 'visible', 'off');
        set(l(3), 'visible', 'off');

        qq_sigma = [qq_sigma, l(1)];
    end

    % combined methods
    qq_sigma_comb = [];
    qq_rms_sigma_comb = nans(1, num_combined_methods);
    for method_index = 1:num_combined_methods
        y = abs(err_est_comb{method_index}(1:1000, :))';
        % y = err_est_comb{method_index}';

        % make plot
        l = qqplot(x(1:num_skip:end), y(1:num_skip:end));
        % adjust plot
        % set(l(1), 'marker', 'none'); %symbols{method_index});
        % set(l(1), 'linestyle', line_symbols{method_index});
        % set(l(1), 'color', colors(2,:));

        set(l(1), 'marker', symbols{method_index});
        set(l(1), 'markeredgecolor', colors(2, :));
        set(l(1), 'linewidth', 0.2);
        % set(l(1), 'markersize', marker_size);
        set(l(1), 'markersize', marker_size);

        set(l(2), 'visible', 'off');
        set(l(3), 'visible', 'off');    

        qq_sigma_comb = [qq_sigma_comb, l(1)];
    end

    box off
    axis equal
    % axis([0 max_error_threshold 0 max_error_threshold])
    axis([0 0.1 0 0.1])
    xlabel('True Error (pix.)', 'fontsize', 16)
    ylabel('Estimated Error (pix.)', 'fontsize', 16)
    % title('Quantile-Quantile Plot', 'fontsize', 16)
    % legend([qq_sigma, qq_sigma_comb], ...
    %     {' IM', ' MC', ' CS', ' Unwt', ' Covar', ' Entropy'}, 'location', 'northoutside', 'NumColumns', 2)
    legend([qq_sigma, qq_sigma_comb], ...
        {individual_method_names{:}, combined_method_names{:}}, 'location', 'northoutside', 'NumColumns', 2)

    % adjust figure size
    set(gcf, 'resize', 'off');
    drawnow();
    set(gcf, 'units', 'inches', 'Position', [500   423   618   518]/user_screen_resolution);
    % set(gca, 'units', 'pix', 'fontsize', 11);
    drawnow();
end
function [qq_rms_indiv, qq_rms_comb] = qq_plot_err_true_est(err, err_est_indv, err_est_comb, individual_method_names, combined_method_names, ...
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
    x = err';
    % only retain finite error values
    x = x(isfinite(x));
    % array to hold plot lines
    l_all = [];
    % plot marker size
    marker_size = 4;
    qq_rms_indiv = nans(1, num_individual_methods);
    qq_rms_comb = nans(1, num_combined_methods);
    % loop through method types (individual and combined)
    for method_type = 1:2
        % extract individual method data
        if method_type == 1
            % extract errors
            err_est = err_est_indv;
            % number of methods
            num_methods = num_individual_methods;
        % extract combined method data
        else
            % extract errors
            err_est = err_est_comb;
            % number of methods
            num_methods = num_combined_methods;
        end

        % loop through method names
        for method_index = 1:num_methods
            % extract estimated error
            % y = err_est{method_index}';
            y = err_est(:, method_index);
            
            % plot
            % l = qqplot(x(1:num_skip:end), y(1:num_skip:end));
            l = qqplot(x, y, 1:100);
            % l = qqplot(abs(x(1:num_skip:end)), abs(y(1:num_skip:end)));

            % calculate rms deviation
            if method_type == 1
                qq_rms_indiv(method_index) = rms(l(1).YData - l(1).XData, 'omitnan');
            else                
                qq_rms_comb(method_index) = rms(l(1).YData - l(1).XData, 'omitnan');
            end

            % adjust plot
            set(l(1), 'marker', symbols{method_index});
            set(l(1), 'markeredgecolor', colors(method_type, :));
            set(l(1), 'linewidth', 0.2);
            set(l(1), 'markersize', marker_size);
            
            set(l(2), 'visible', 'off');
            set(l(3), 'visible', 'off');
            
            % aggregate handles to lines
            l_all = [l_all, l(1)];                        
        end
    end

    % adjust plot
    box off
    axis equal    
    % axis([-0.1 0.1 -0.1 0.1])
    axis([-max_error_threshold, max_error_threshold, -max_error_threshold, max_error_threshold])
    set(gca, 'YAxisLocation', 'origin')
    set(gca, 'XAxisLocation', 'origin')
    % move_axes_to_origin(gca);
    % axis([0 0.1 0 0.1])
    set(gca, 'xtick', -0.1:0.1:0.1)
    % xl = xlabel('True Error (pix.)', 'fontsize', 10); %16)
    xl = xlabel('\epsilon_{True} (pix.)', 'fontsize', 10); %16)    
    set(xl, 'HorizontalAlignment', 'center');

    % yl = ylabel('Estimated Error (pix.)', 'fontsize', 10); %16)
    yl = ylabel('\epsilon_{Est} (pix.)', 'fontsize', 10); %16)
    yl.Position(1) = -0.02;
    set(yl, 'HorizontalAlignment', 'right');
    legend(l_all, ...
        {individual_method_names{:}, combined_method_names{:}}, 'location', 'northoutside', 'orientation', 'horizontal') %, 'NumColumns', 2)
    % adjust figure size
    set(gcf, 'resize', 'off');
    drawnow();
    % set(gcf, 'units', 'inches', 'Position', [500   423   618   518]/user_screen_resolution);
    set(gcf, 'units', 'inches', 'Position', [500  541  400  400]/user_screen_resolution);
    
    drawnow();

end
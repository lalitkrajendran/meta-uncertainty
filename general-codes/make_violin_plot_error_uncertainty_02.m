function violins = make_violin_plot_error_uncertainty_02(err_abs, sigma_indv, sigma_comb, err_rms, sigma_rms_indv, sigma_rms_comb, ...
                                                            bins, individual_method_array, combination_method_array, ...
                                                            colors, user_screen_resolution, max_error_threshold)

% Function to make violin plots of error and uncertainty
%
% INPUTS: (all in pixels)
% err_abs: absolute error
% sigma_indv, sigma_comb: individual and combined uncertainties
% err_rms: rms of error
% sigma_rms_indv, sigma_rms_comb: rms of individual and combined uncertainties
% bins: histogram bins
% individual_method_array, combination_method_array: names of individual and combined uncertainty methods
% colors: colors for individual and combined uncertainty schemes
% user_screen_resolution: dpi
% max_error_threshold: max allowable error
%
% OUTPUTS:
% NONE
%
% AUTHOR:
% Lalit Rajendran (lrajendr@purdue.edu)

    % ==========================
    % number of methods
    % ==========================
    num_individual_methods = numel(individual_method_array);
    num_combined_methods = numel(combination_method_array);
    num_methods_total = num_individual_methods + num_combined_methods + 1;

    % ensure absoluve value of error
    err_abs = abs(err_abs);

    % --------------------------
    % obtain max number of elements in any of the arrays
    % --------------------------
    N_all = nans(1, num_individual_methods + num_combined_methods + 1);    
    N_all(1) = numel(err_abs);
    for i = 1:num_individual_methods
        N_all(1 + i) = numel(sigma_indv{i});
    end
    for i = 1:num_combined_methods
        N_all(1 + num_individual_methods + i) = numel(sigma_comb{i});
    end
    N_max = max(N_all);

    % --------------------------
    % pad arrays with nans to ensure they are the same size
    % --------------------------
    err_abs = padarray(err_abs, [0, N_max - N_all(1)], NaN, 'post');
    for method_index = 1:num_individual_methods
        sigma_indv{method_index} = padarray(sigma_indv{method_index}, [0, N_max - N_all(method_index+1)], NaN, 'post');
    end
    for method_index = 1:num_combined_methods
        sigma_comb{method_index} = padarray(sigma_comb{method_index}, [0, N_max - N_all(method_index+1+num_individual_methods)], NaN, 'post');
    end

    % --------------------------
    % make violin plot arrays and colors
    % --------------------------
    plot_array = nans(N_max, num_individual_methods + num_combined_methods + 1);
    colors_array = nans(num_methods_total, 3);
    for i = 1:num_individual_methods
        plot_array(:, i) = sigma_indv{i}';
        colors_array(i, :) = colors(1, :);
    end

    plot_array(:, num_individual_methods+1) = err_abs';
    colors_array(num_individual_methods+1, :) = [0, 0, 0];

    for i = 1:num_combined_methods
        plot_array(:, i + num_individual_methods + 1) = sigma_comb{i}';
        colors_array(i + num_individual_methods + 1, :) = colors(2, :);
    end

    % --------------------------
    % plot violins
    % --------------------------
    h = figure;
    violins = violin_plot_pdf(plot_array, bins, {individual_method_array{:}, 'Error', combination_method_array{:}}, colors_array);
    % drawnow();
    figure(h)
    % face colors
    color_all = cell(1, numel(violins));

    % ==========================
    % adjust violin colors and add lines corresponding to rms
    % ==========================
    for violin_index = 1:numel(violins)
        % x co-ordinates of current violin plot
        x = violins(violin_index).XData;

        % --------------------------
        % individual uncertainty methods
        % --------------------------
        if violin_index <= num_individual_methods
            % adjust color
            violins(violin_index).FaceColor = colors(1, :);
            drawnow();
            % extract rms
            y = sigma_rms_indv(violin_index);
            
        % --------------------------
        % error
        % --------------------------
        elseif violin_index == num_individual_methods+1
            % adjust color
            violins(violin_index).FaceColor = [0, 0, 0];
            drawnow();
            % extract rms
            y = err_rms;
            plot([min(violins(1).XData), max(violins(numel(violins)).XData)], y*[1, 1], '--', 'color', [0, 0, 0, 0.2])
            drawnow();
        
        % --------------------------
        % combined uncertainty methods
        % --------------------------
        elseif violin_index > num_individual_methods + 1
            % adjust color
            violins(violin_index).FaceColor = colors(2, :);
            drawnow();
            % extract rms            
            y = sigma_rms_comb(violin_index - (num_individual_methods + 1));
        end

        % --------------------------
        % adjust violin properties
        % --------------------------
        color_all{violin_index} = violins(violin_index).FaceColor;
        violins(violin_index).FaceAlpha = 0.25;
        drawnow();

        % --------------------------
        % plot rms value
        % --------------------------
        % plot([min(x), max(x)], y*[1, 1], 'color', [violins(violin_index).ViolinColor, violins(violin_index).ViolinAlpha+0.2], 'linewidth', 3)
        plot([min(x), max(x)], y*[1, 1], 'color', [violins(violin_index).FaceColor, violins(violin_index).FaceAlpha+0.2], 'linewidth', 3)
        drawnow();

        % --------------------------
        % add text for rms value
        % --------------------------
        x_text = mean(x) + 0.1 * (max(x) - min(x));        
        y_text = y + 0.01;
        text(x_text, y_text, num2str(y, '%.2f'), 'fontweight', 'bold', 'color', violins(violin_index).FaceColor); 
        
    end

    % annotate y axis
    ylabel('Error and Uncertainty (pix.)', 'fontsize', 14)
    % adjust limits
    ylim([0 max_error_threshold])
    drawnow();
    
    % turn off x axis line
    ax = gca(h);
    ax.XAxis.Axle.Visible = 'off';
    ax.XAxis.TickLength = [0 0];
    
    % turn off y axis line
    ax.YAxis.Axle.Visible = 'off';
    ax.YAxis.TickLength = [0.005 0.005];

    % adjust figure position
    set(h, 'resize', 'off');
    drawnow();
    set(h, 'units', 'inches', 'Position', [352   526   895   392]/user_screen_resolution)
    drawnow();
    
    % set(gca, 'units', 'pix', 'fontsize', 11);    
end

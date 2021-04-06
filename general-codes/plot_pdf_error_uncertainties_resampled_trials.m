function plot_pdf_error_uncertainties_resampled_trials(pdf_err_resampled_trials, pdf_unc_resampled_trials, bins, method_names, skip_trials, user_screen_resolution)
% Function to plot pdf of errors and uncertainties for each set of resampling trials
%
% INPUTS:
% pdf_err_resampled_trials: pdf of errors from resampling
% pdf_unc_resampled_trials: pdf of uncertainties from resampling
% bins: histogram bins
% method_names: name of the uncertainty schemes
% skip_trials: number of trials to skip
% user_screen_resolution: desired dpi
%
% OUTPUTS:
% None
%
% AUTHOR:
% Lalit Rajendran (lrajendr@purdue.edu)
%
% DATE:
% 2020/06/19

    % number of methods
    num_methods = numel(method_names);
    % get line colors
    colors = lines(1);
    % number of valid trials
    num_valid_trials = size(pdf_unc_resampled_trials{1}, 1);
    
    figure
    for method_index = 1:num_methods+1
        subplot(num_methods+1, 1, method_index)
        if method_index == 1
            pdf_current = pdf_err_resampled_trials;
        else
            pdf_current = pdf_unc_resampled_trials{method_index-1};
        end

        for trial_index = 1:skip_trials:num_valid_trials
            plot(bins(1:end-1), pdf_current(trial_index, :), 'color', [colors(1, :) 0.5])
            hold on
        end
        xlim([0 bins(end)])
        box off
        if method_index == 1
            title('Error')
        else
            title(method_names{method_index-1})
        end
        if method_index == num_methods+1
            xlabel('(pix.)')
        else
            set(gca, 'xticklabel', []);
        end
    end

    % set common y limit
    % set_common_ylim(gcf, 'manual', [0, max(y_max)]);
    set_common_ylim(gcf, 'auto');

    % adjust plot
    set(gcf, 'resize', 'off');
    % set(gcf, 'position', [680   370   530   480]);
    set(gcf, 'units', 'inches', 'position', [680   370   530   480]/user_screen_resolution);
    set(gca, 'units', 'pix', 'fontsize', 11);

end
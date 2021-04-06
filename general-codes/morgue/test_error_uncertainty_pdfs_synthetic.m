clear
close all
clc

addpath(genpath('/scratch/shannon/c/aether/Projects/BOS/general-codes/matlab-codes/'));
setup_default_settings;
%%
sigma_mean = 0.1;
sigma_std = 0.05;

num_trials = 1e6;
max_error_threshold = 0.5;

bins = linspace(0, max_error_threshold, sqrt(num_trials));

sigma_all = nans(1, num_trials);
err_all = nans(1, num_trials);


parfor trial_index = 1:num_trials
    % random uncertainty
    sigma_all(trial_index) = sigma_mean + randn() * (sigma_std - sigma_mean);
    if sigma_all(trial_index) <= 0 
	sigma_all(trial_index) = NaN;
	continue;
    end
    % randomly sample error
    err_all(trial_index) = randn() * sigma_all(trial_index);
end

%%
sigma_all_valid = sigma_all(sigma_all > 0);
err_all_valid = err_all(abs(err_all) < max_error_threshold);

%% pdf
[N_err, edges_err] = histcounts(abs(err_all), bins, 'Normalization', 'pdf');
[N_sigma, edges_sigma] = histcounts(abs(sigma_all), bins, 'Normalization', 'pdf');

cdf_err = histcounts(abs(err_all), bins, 'Normalization', 'cdf');
cdf_sigma = histcounts(sigma_all, bins, 'Normalization', 'cdf');

rms_err = rms(err_all(isfinite(err_all)));
rms_sigma = rms(sigma_all(isfinite(sigma_all)));

%% plot
colors = lines(2);
y_max = max([N_err, N_sigma]);
figure
subplot(2, 1, 1)
plot(edges_err(1:end-1), N_err, 'color', colors(1, :))
hold on
plot([rms_err, rms_err], [0, y_max], '--', 'color', colors(1, :))

plot(edges_sigma(1:end-1), N_sigma, 'color', colors(2, :))
plot([rms_sigma, rms_sigma], [0, y_max], '--', 'color', colors(2, :))

xlim([0, max_error_threshold])
box off
set(gca, 'xticklabel', '')
title('PDF')
legend('Error', 'RMS Error', 'Uncertainty', 'RMS Uncertainty')

subplot(2, 1, 2)
plot(edges_err(1:end-1), cdf_err, 'color', colors(1, :))
hold on
plot(edges_sigma(1:end-1), cdf_sigma, 'color', colors(2, :))

xlim([0, max_error_threshold])
box off
xlabel('(pix.)')
title('CDF')
legend('Error', 'Uncertainty', 'location', 'southeast')

set(gcf, 'Position', [680   751   457   648])

%% qqplot of error vs uncertainty

line_symbols = {'-'; ':'; '-.'};
figure
plot([0 max_error_threshold], [0 max_error_threshold], 'k');
hold on

% make plot
l = qqplot(err_all, sigma_all);
% adjust plot
%     set(l(1), 'marker', 'none');
%     set(l(1), 'linestyle', line_symbols{method_index});
%     set(l(1), 'color', colors(1,:));

set(l(1), 'marker', 'o');
set(l(1), 'markeredgecolor', colors(1, :));
set(l(1), 'markersize', 4);

set(l(2), 'visible', 'off');
set(l(3), 'visible', 'off');

box off
axis equal
axis([0 max_error_threshold 0 max_error_threshold])
xlabel('Error (pix.)', 'fontsize', 16)
ylabel('Uncertainty (pix.)', 'fontsize', 16)
title('Quantile-Quantile Plot', 'fontsize', 16)

set(gcf, 'resize', 'off');
pause(0.1);
set(gcf, 'Position', [716   423   618   518]);

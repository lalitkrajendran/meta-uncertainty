clear
close all
clc

restoredefaultpath;
addpath(genpath('/scratch/shannon/c/aether/Projects/BOS/general-codes/matlab-codes/'));
addpath('../prana-uncertainty-average-dc-new-im-2-save-deform-cs/');
addpath('../stereo_uncertainty_codes_packaged/');
setup_default_settings;

N = 1e3;

x = randn(1, N);
y = randn(1, N);
% y = x;

num_bins = sqrt(N);

[cdfxy, edges_x, edges_y] = histcounts2(x, y, num_bins, 'normalization', 'cdf');

% q_levels = [0.25, 0.75];

M = contourc(edges_x(1:end-1), edges_y(1:end-1), cdfxy, 0:0.05:1);

c = find(M(1, :) == 0.25);
num_c = M(2, c);
xarr1 = M(1, (c+1):(c+1+num_c-1));
yarr1 = M(2, (c+1):(c+1+num_c-1));

c = find(M(1, :) == 0.75);
num_c = M(2, c);
xarr3 = M(1, (c+1):(c+1+num_c-1));
yarr3 = M(2, (c+1):(c+1+num_c-1));

x_all = [xarr1, fliplr(xarr3), xarr1(1)];
y_all = [yarr1, fliplr(yarr3), yarr1(1)];

[iqr_bivar, ind1, ind2] = calculate_minimum_distance_between_two_curves(xarr1, yarr1, xarr3, yarr3);
iqr_bivar = polyarea(x_all, y_all) - iqr(x) * iqr(y);

fprintf('iqr: %.2f\n', iqr_bivar);

figure
subplot(1, 2, 1)
surf(edges_x(1:end-1), edges_y(1:end-1), cdfxy, 'linestyle', 'none')
colormap('gray')
hold on
plot3(xarr1, yarr1, 0.25 * ones(size(xarr1)), 'r')  
plot3(xarr3, yarr3, 0.75 * ones(size(xarr3)), 'b')
xlim([edges_x(1), edges_x(end)])
ylim([edges_y(1), edges_y(end)])
zlim([0 1])
xlabel('x')
ylabel('y')
zlabel('CDF')

subplot(1, 2, 2)
plot(xarr1, yarr1, 'ro-', xarr3, yarr3, 'bo-', 'linewidth', 2);
hold on
plot([xarr1(ind1), xarr3(ind2)], [yarr1(ind1), yarr3(ind2)], 'go-')
fill(x_all, y_all, 'k', 'facealpha', 0.2, 'edgecolor', 'none')
axis equal
xlim([min(x_all), max(x_all)])
ylim([min(y_all), max(y_all)])
box off
xlabel('x')
ylabel('y')

set(gcf, 'resize', 'off')
set(gcf, 'position', [142         256        1125         499])

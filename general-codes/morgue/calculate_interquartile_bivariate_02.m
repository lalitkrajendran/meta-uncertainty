% Code to test quartile calculation for bivariate data

clear
close all
clc

restoredefaultpath;
addpath(genpath('/scratch/shannon/c/aether/Projects/BOS/general-codes/matlab-codes/'));
addpath('../prana-uncertainty-average-dc-new-im-2-save-deform-cs/');
addpath('../stereo_uncertainty_codes_packaged/');
addpath ../CompPD
addpath ../MvLogNRand
setup_default_settings;

% ======================
% random variable settings
% ======================
% number of elements in the random vector
N = 1e3;
% mean
mu = [0, 0];
% mu = [11, 12];
% standard deviation
sigma = [1 0.3; 0.3, 1];
% sigma = [.1 .3];
% correlation coefficient matrix
% CorrMat = [1 .2 .4 ; .2 1 .5 ; .4  .5 1];
% CorrMat = [1 .2 ; .2 1];
% desired depth value for contours
alpha0 = 0.5;

% ======================
% generate random variable
% ======================
% create bivariate random variable
X = mvnrnd(mu, sigma, N);
% X = MvLogNRand( mu , sigma , N , CorrMat );
% extract size of X
[n, p] = size(X);

% ======================
% calculate projection depth
% ======================
% Obtain the optimal direction vectors for approximate computing the projection depth 
% and its associated estimators
AppVec1 = AppVecPD(X, 1e2);
% calculate depth values
pdv = PDVal(X, AppVec1);

% ======================
% calculate depth weighted statistics
% ======================
C = median(pdv);
K = 3;
% calculate median
pm0 = PM(X, AppVec1, true);
% calculate projection weighted standard deviation
[pws, pm0, weitss, tmpMat] = PWS(X, pdv, C, K, C, K);
fprintf('pws: %.2f\n', pws);
% calculate ordinary covariance matrix
covx = cov(X);
% calculate contours of constant depth values
vpmat= PC2D(X, AppVec1, alpha0, true, false);
% calculate area inside the contour region
iqr_area = polyarea(vpmat{1}(:, 1), vpmat{1}(:, 2));
fprintf('iqr_area: %.2f\n', iqr_area);

% ======================
% plot results
% ======================

% Plot the scatter plot of X
figure(1)
subplot(1, 2, 1)
plot(X(:, 1), X(:, 2), '.')
xlabel('X_{1}'); ylabel('X_{2}');
box on;
hold on
% plot median
plot(pm0(1), pm0(2), 'r+');
% plot depth contours
for i = 1:numel(vpmat)
    plot(vpmat{i}(:, 1), vpmat{i}(:, 2), 'g')
end
axis equal

colors = cbrewer('seq', 'Blues', 100);
subplot(1, 2, 2)
hold on
for i = 1:N
    color_index = weitss(i)/max(weitss) * 99;
    % plot(X(i, 1), X(i, 2), '.', 'color', [colors(1, :), weitss(i)/max(weitss)])
    plot(X(i, 1), X(i, 2), '.', 'color', colors(round(color_index) + 1, :))
end
% plot median
plot(pm0(1), pm0(2), 'r+');
% plot depth contours
for i = 1:numel(vpmat)
    plot(vpmat{i}(:, 1), vpmat{i}(:, 2), 'g')
end
axis equal
box on;
xlabel('X_{1}'); ylabel('X_{2}');
title(sprintf('C: %.2f, K: %.2f', C, K));
set(gcf, 'resize', 'off');
set(gcf, 'position', [200   334   922   426]);
set(gcf, 'resize', 'off');
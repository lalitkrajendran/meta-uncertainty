clear
close all
clc

restoredefaultpath;
addpath(genpath('/scratch/shannon/c/aether/Projects/BOS/general-codes/matlab-codes/'));
setup_default_settings;

% ===========
% settings
% ===========
% number of elements
N = 1e6;
% mean
mu = [0 0];
% standard deviation
% sigma = [1 1.5; 1.5 3];
sigma = [1 0.5; 0.5 3];
% bins
% bins = linspace(-10, 10, sqrt(N));
bins = linspace(-10, 10, 100);

% generate random numbers
X = mvnrnd(mu, sigma.^2, N);
% X = randn(N, 2);

% ===========
% calculate probabilities
% ===========
% individual
p1 = histcounts(X(:, 1), bins, 'normalization', 'probability');
p2 = histcounts(X(:, 2), bins, 'normalization', 'probability');

pdf1 = histcounts(X(:, 1), bins, 'normalization', 'pdf');
pdf2 = histcounts(X(:, 2), bins, 'normalization', 'pdf');

% joint
p12 = histcounts2(X(:, 1), X(:, 2), 'xbinedges', bins, 'ybinedges', bins, 'normalization', 'probability');

pdf11 = histcounts2(X(:, 1), X(:, 1), 'xbinedges', bins, 'ybinedges', bins, 'normalization', 'pdf');
pdf22 = histcounts2(X(:, 2), X(:, 2), 'xbinedges', bins, 'ybinedges', bins, 'normalization', 'pdf');
pdf12 = histcounts2(X(:, 1), X(:, 2), 'xbinedges', bins, 'ybinedges', bins, 'normalization', 'pdf');

% ===========
% calculate entropies from probabilitles
% ===========
% calculate individual entropies
e1 = -sum(p1 .* log2(p1), 'omitnan');
e2 = -sum(p2 .* log2(p2), 'omitnan');

% calculate joint entropies 
e12 = -sum(p12 .* log2(p12), 'all', 'omitnan');

% calculate mutual information directly
m12_d = -sum(p12 .* log2(p12./p1 * p2'), 'all', 'omitnan');
% calculate mutual information
m12 = e1 + e2 - e12;

% ===========
% calculate entropies from pdfs
% ===========
% calculate individual entropies 
h = pdf1 .* log(pdf1);
h(isnan(h)) = 0;
e_pdf_1 = -trapz(bins(1:end-1), h);
e_pdf_1b = calculate_self_entropy(X(:, 1), bins);

h = pdf2 .* log(pdf2);
h(isnan(h)) = 0;
e_pdf_2 = -trapz(bins(1:end-1), h);
e_pdf_2b = calculate_self_entropy(X(:, 2), bins);

% joint
h = pdf12 .* log(pdf12);
h(isnan(h)) = 0;
e_pdf_12 = -trapz(bins(1:end-1), trapz(bins(1:end-1), h));
e_pdf_12b = calculate_joint_entropy(X(:, 1), X(:, 2), bins, bins);

% calculate mutual information
m_pdf_12 = e_pdf_1 + e_pdf_2 - e_pdf_12;
m_pdf_12b = calculate_mutual_information(X(:, 1), X(:, 2), bins, bins);

% calculate theoretical mutual information 
m_pdf_12_theory = -1/2 * log(1 - sigma(1, 2)^2);

% ===========
% display results
% ===========
fprintf('sigma: %.2f, %.2f, %.2f\n', sigma(1, 1), sigma(2, 2), sigma(1, 2));
fprintf('entropy: %.2f, %.2f, %.2f\n', e1, e2, e12);
fprintf('entropy, pdf: %.2f, %.2f, %.2f\n', e_pdf_1, e_pdf_2, e_pdf_12);
fprintf('entropy, pdf, b: %.2f, %.2f, %.2f\n', e_pdf_1b, e_pdf_2b, e_pdf_12b);
fprintf('mutual information: %2f, %.2f, %.2f\n', m12, m_pdf_12, m_pdf_12b);



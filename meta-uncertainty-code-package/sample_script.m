% This script performs resampling calculations for a range of particle removal
% percentage and plots the change in uncertainties for different methods

clear
close all
clc

restoredefaultpath;
addpath(genpath('/scratch/shannon/c/aether/Projects/BOS/general-codes/matlab-codes/'));
addpath(genpath('/scratch/shannon/c/aether/Projects/BOS/general-codes/piv-image-generation/'));
addpath('../prana/')
addpath('../general-codes/')
addpath('../histogram_distance/')

setup_default_settings;

% ==========================
% settings
% ==========================
% processing pass number
pass_number = 4;
% array of uncertainty methods
individual_method_array = {'IM'; 'MC'; 'CS'};
num_individual_methods = numel(individual_method_array);
% array of snr metrics
snr_metric_array = {'PPR'; 'MI'};
num_snr_methods = numel(snr_metric_array);
% components
components = {'x'; 'y'};
num_components = numel(components);

% --------------------------
%% resampling settings
% --------------------------
% resampling method names
resampling_method_name_array = {'remove-paired'; 'remove-random'; 'add-random'};
% resampling_method_name_array = {'add-random'};
num_resampling_methods = numel(resampling_method_name_array);
% number of particles to remove
percentage_particles_remove_array = 0.05:0.05:0.25; %0.025:0.025:0.25; 
num_ppr = numel(percentage_particles_remove_array);
% number of resampling trials
num_resampling_trials = 1e1;
% display_resampled_images? (true/false)
display_resampled_images = 0;

% --------------------------
%% particle identification settings
% --------------------------
% intensity threshold for particle identification
intensity_threshold_percentile = 10;
% particle diameter (pix.)
d_p = 3;
% particle sizing settings
sizeprops = struct;
sizeprops.method = 'IWC';
sizeprops.p_area = 2; % this will be modified for each image
sizeprops.sigma = 4;
sizeprops.errors = 0;
sizeprops.save_dir = [];
% display id results?
display_id_results = 0;

% --------------------------
%% weight settings
% --------------------------
% number of bins for weights
num_bins_w = 30;
% bins for weights
bins_w = linspace(0, 1, num_bins_w);
% metric type ('rms', 'iqr', 'std')
metric_name = 'iqr';

% ==========================
% load data
% ==========================
% load jobfile
jobfile = load('jobfile.mat');
jobfile = jobfile.Data;

% load vector field
results = load('sample-displacement.mat');

% load deformed images
saved_img = load('im1-def.mat');
im1d = saved_img.im1d;
saved_img = load('im2-def.mat');
im2d = saved_img.im2d;

% select a grid point row and column
r = 10; 
c = 10;

% ================================================
%% allocate memory
% ================================================
fprintf('allocating memory for variables\n');

% uncertainties
unc_resampling = cell(num_ppr, num_resampling_methods);

% snr metric
snr_resampling = cell(num_ppr, num_resampling_methods);                        

% ==========================
% run resampling
% ==========================
fprintf('running monte-carlo\n');

% set seed for random number generator
rng();

% start timer
tic

% --------------------------
% extract job file properties
% --------------------------
% extract pass settings
current_pass_settings = jobfile.(['PIV' num2str(pass_number)]);
% extract window size
[window_size_x, window_size_y] = extract_window_size(current_pass_settings);
% extract window resolution
[window_resolution_x, window_resolution_y] = extract_window_resolution(current_pass_settings);        
% extract grid resolution
[grid_resolution_x, grid_resolution_y] = extract_grid_resolution(current_pass_settings);

% --------------------------
%% extract processing results for this grid point
% --------------------------
% extract co-ordinates
X = results.X(r, c);
Y = results.Y(r, c);

% extract individual uncertainties for this grid point
unc_base = extract_planar_uncertainties(results.uncertainty2D, r, c);

% extract snr metric
snr_base = extract_snr_metric(results.SNRmetric, r, c);

% --------------------------
% calculate base uncertainty
% --------------------------    
% extract image height and width
[image_height, image_width] = size(im1d);

% extract interrogation window
im1_sub = extract_interrogation_window(im1d, X, Y, [window_size_x, window_size_y]);
im2_sub = extract_interrogation_window(im2d, X, Y, [window_size_x, window_size_y]);

% calculate uncertainties of this interrogation window
[unc_sub, snr_sub] = calculate_planar_uncertainties(im1_sub, im2_sub, current_pass_settings, individual_method_array);

% calculate average particle diameter
d_p = results.uncertainty2D.Autod(r, c)/sqrt(2);

% --------------------------
% loop through particle removal percentages
% --------------------------
for particle_remove_index = 1:num_ppr
    % current percentage to remove
    percentage_particles_remove = percentage_particles_remove_array(particle_remove_index);
    fprintf('particle remove percentage: %.2f\n', percentage_particles_remove);
    
    % --------------------------
    % loop through resampling methods
    % --------------------------
    for resampling_method_index = 1:num_resampling_methods
        fprintf('resampling method index: %d\n', resampling_method_index);
        resampling_method_name = resampling_method_name_array{resampling_method_index};
        % calculate uncertainty from resampling
        [unc, snr] = calculate_uncertainty_resampling_general(im1_sub, im2_sub, d_p, sizeprops, resampling_method_name, percentage_particles_remove, ...
                                                                num_resampling_trials, individual_method_array, intensity_threshold_percentile, ...
                                                                current_pass_settings, display_resampled_images);
        
        % copy results
        unc_resampling{particle_remove_index, resampling_method_index} = unc;
        snr_resampling{particle_remove_index, resampling_method_index} = snr;
    end    
end

% ==========================
% calculate weights
% ==========================
% --------------------------
% initialize 
% --------------------------            
unc_ratio = cell(num_resampling_methods, 1);
unc_ratio_fit = cell(1, num_resampling_methods, 1);
unc_ratio_rate = cell(1, num_resampling_methods, 1);
weights = cell(num_resampling_methods, 1);        

% --------------------------
% weight calculation
% --------------------------
for resampling_method_index = 1:num_resampling_methods              
    [weights{resampling_method_index}, unc_ratio{resampling_method_index}, unc_ratio_fit{resampling_method_index}] = calculate_weights_from_resampled_uncertainty_ratios(unc_sub, {unc_resampling{:, resampling_method_index}}, ...
                                                    metric_name, individual_method_array, percentage_particles_remove_array, components, num_resampling_trials, 'struct');
end    

% ==========================
% calculate combined uncertainty
% ==========================
% --------------------------
% extract individual uncertainty
% --------------------------
unc_indiv.x = nans(1, num_individual_methods);
unc_indiv.y = nans(1, num_individual_methods);

for individual_method_index = 1:num_individual_methods
    method_name = lower(individual_method_array{individual_method_index});
    unc_indiv.x(individual_method_index) = unc_base.([method_name 'x']);
    unc_indiv.y(individual_method_index) = unc_base.([method_name 'y']);
end

% --------------------------
% calculate combined uncertainty
% --------------------------
unc_comb.x = nans(1, num_resampling_methods);
unc_comb.y = nans(1, num_resampling_methods);
for resampling_method_index = 1:num_resampling_methods
    unc_comb.x(resampling_method_index) = sum(unc_indiv.x .* weights{resampling_method_index}.x, 2); 
    unc_comb.y(resampling_method_index) = sum(unc_indiv.y .* weights{resampling_method_index}.y, 2); 
end                

% shut down parallel pool
shut_down_parpool();

% stop timer
toc

% --------------------------
% save results
% --------------------------
% file name
% filename = ['resampling-statistics-nr='  num2str(num_resampling_trials, '%d') '.mat'];
filename = 'sample-results.mat';
% save results to file
save(filename, 'r', 'c', 'X', 'Y', 'unc_sub', 'snr_sub', 'unc_base', 'snr_base', ...
                'unc_sub', 'snr_sub', 'unc_resampling', 'snr_resampling', ...
                'weights', 'unc_ratio', 'unc_ratio_fit', 'unc_indiv', 'unc_comb');   

% This script calculates the pdf of uncertainties estimated by three
% different algorithms for experimental images by jack-knife resampling.

%%
clear
close all
clc

restoredefaultpath;
addpath(genpath('/scratch/shannon/c/aether/Projects/BOS/general-codes/matlab-codes/'));
addpath('../prana-uncertainty-average-dc-new-im-2-save-deform-cs/');
setup_default_settings;
dbstop if error

%% read/write settings

% window resolution
window_resolution_array = [32, 64];
num_window_size = numel(window_resolution_array);

% pass number
pass_number = 4;

% directory containing images
top_image_directory = fullfile('/scratch/shannon/a/aether/Projects/PIV_PTV_Uncertainty_Quantification/Planar_uncertainty_work/Images/');
% directory containing files to be read
top_read_directory = fullfile('/scratch/shannon/a/aether/Projects/PIV_PTV_Uncertainty_Quantification/uncertainty-combination/analysis/results/', 'experiment-new');

% array of data set names
dataset_name_array = {'PivChal03B'; 'PivChal05B'; 'stagnation_flow'; 'Vortex_Ring'; 'Jetdata'};
num_datasets = numel(dataset_name_array);

% array of uncertainty methods
uncertainty_method_array = {'IM'; 'MC'; 'CS'};
num_uncertainty_methods = numel(uncertainty_method_array);

% directory where results of this analysis are to be saved
top_write_directory = fullfile('/scratch/shannon/a/aether/Projects/PIV_PTV_Uncertainty_Quantification/uncertainty-combination/analysis/results/', 'monte-carlo', 'new-processing-all-datasets');
mkdir_c(top_write_directory);

% number of trials
num_trials = 1e4;

%% resampling settings

% number of particles to remove
% num_particles_remove_array = [5, 10, 15];
percentage_particles_remove = 0.25;
% number of resampling trials
num_resampling_trials = 1e3;

% uncertainty flags for the various methods (0 = False, 1 = True)
uncertainty_flags.ppruncertainty = 0;
uncertainty_flags.miuncertainty = 0;
uncertainty_flags.mcuncertainty = 1;

%% particle identification settings

% intensity threshold for particle identification
intensity_threshold = 10;

% particle diameter (pix.)
d_p = 3;

% particle sizing settings
sizeprops = struct;
sizeprops.method = 'IWC';
sizeprops.p_area = 2; % this will be modified for each image
sizeprops.sigma = 4;
sizeprops.errors = 0;
sizeprops.save_dir = [];

%% plot settings

% save_figure? (true/false)
save_figures = true;

%% directory settings for this case

% directory to save results for this case
current_write_directory = fullfile(top_write_directory, ['trials=' num2str(num_trials, '%d')], ['ppr' num2str(percentage_particles_remove, '%.2f') '_nrs' num2str(num_resampling_trials)], 'new');
mkdir_c(current_write_directory);

% directory to save figures
figure_write_directory = fullfile(current_write_directory, 'figures');
mkdir_c(figure_write_directory);

%% pre-load all dataset

fprintf('Loading all datasets into memory\n');

results_all = cell(num_window_size, num_datasets);
jobfile_all = cell(num_window_size, num_datasets);
files_im1 = cell(num_window_size, num_datasets);
files_im2 = cell(num_window_size, num_datasets);

% loop through window resolutions
for window_size_index = 1:num_window_size
    % loop through datasets
    for dataset_index = 1:numel(dataset_name_array)
        % name of the data set
        dataset_name = dataset_name_array{dataset_index};
        
        fprintf('Dataset: %s\n', dataset_name);

        %% Load data

        % results directory for current data set
        current_results_directory = fullfile(top_read_directory, dataset_name, ['WS' num2str(window_size_index)]);

        % directory containing vectors
        vectors_directory = fullfile(current_results_directory, 'vectors');

        % load results for vectors, errors and uncertainties
        results_all{window_size_index, dataset_index} = load_directory_data(vectors_directory, ['*pass' num2str(pass_number) '*.mat']);

        % load jobfile
        jobfile_all{window_size_index, dataset_index} = load(fullfile(current_results_directory, 'jobfile.mat'));
        
        %% load listing of deformed images
        
        % directory containing deformed images for current data set
        deformed_images_directory = fullfile(top_read_directory, dataset_name, ['WS' num2str(window_size_index)], 'vectors', 'imDeform');
        
        % get list of im1 files in the directory
        files_im1{window_size_index, dataset_index} = get_directory_listing(deformed_images_directory, 'PIV*im1d*.tif');
        % get list of im2 files in the directory
        files_im2{window_size_index, dataset_index} = get_directory_listing(deformed_images_directory, 'PIV*im2d*.tif');
        
    end
end

%% load errors for all datasets

fprintf('Loading all errors into memory\n');

errors_all = cell(num_window_size, num_datasets);

% loop through window resolutions
for window_size_index = 1:num_window_size
    % loop through datasets
    for dataset_index = 1:numel(dataset_name_array)
        % name of the data set        
        dataset_name = dataset_name_array{dataset_index};
        
        fprintf('Dataset: %s\n', dataset_name);
        
        %% Load data

        % results directory for current data set
        current_results_directory = fullfile(top_read_directory, dataset_name, ['WS' num2str(window_size_index)]);

        % load results for vectors, errors and uncertainties
        errors_all{window_size_index, dataset_index} = load_directory_data(current_results_directory, 'errors.mat');
        
    end
end

%% allocate memory

fprintf('allocating memory for variables\n');
% dataset address
dataset_index_array = nans(1, num_trials);
window_resolution_index_array = nans(1, num_trials);

% errors
err_U = nans(1, num_trials);
err_V = nans(1, num_trials);

% uncertainties
unc_trials = cell(1, num_trials);

% uncertainty from resampling
unc_resampling_avg = cell(1, num_trials);
unc_resampling_std = cell(1, num_trials);
unc_resampling = cell(1, num_trials);

%% loop through trials and perform analysis

fprintf('running monte-carlo\n');

% set seed for random number generator
rng(0);

% generate random numbers for overall dataset index
% overall_dataset_index_array = 1 + round(rand(1, num_trials) * (num_datasets * num_window_resolution - 1));
overall_dataset_index_array = 1 + round(rand(1, num_trials) * (num_datasets - 1));

overall_window_resolution_index_array = 1 + round(rand(1, num_trials) * (num_window_size - 1));

parfor trial_index = 1:num_trials
    fprintf('trial_index: %d\n', trial_index);

    % calculate random number for data set
    dataset_index = randi(num_datasets);
    dataset_index_array(trial_index) = dataset_index;
    % calculate random number for window size
    window_resolution_index_array(trial_index) = randi(num_window_size);
    
    % load errors
    errors = errors_all{window_resolution_index_array(trial_index), dataset_index_array(trial_index)};

    % calculate number of snapshots available for the current dataset and
    % window size
    num_snapshots = numel(results_all{window_resolution_index_array(trial_index), dataset_index_array(trial_index)});
    if num_snapshots > size(errors.err_U, 3)
        num_snapshots = size(errors.err_U, 3);
    end
    % calculate random number for the snapshot
    snapshot_index = randi(num_snapshots);

    % load dataset
    results = results_all{window_resolution_index_array(trial_index), dataset_index_array(trial_index)}(snapshot_index);

    % load jobfile
    jobfile = jobfile_all{window_resolution_index_array(trial_index), dataset_index_array(trial_index)};
        
    % calculate size of the results array
    [num_rows, num_cols] = size(results.X);
    
    %% select random number for grid point in this snapshot
    
    if dataset_index <= 3
        c = 1:size(results.X, 2);
        r = 1:size(results.X, 1);
    
    elseif dataset_index == 4
        % extract co-ordinate arrays for the processing grid
        xvec = results.X(1, :);
        yvec = results.Y(:, 1);
        
        % find vector locations in the processing results for which
        % true solution is available
        c = find(xvec > errors.X(1, 1) & xvec < errors.X(1, end));
        r = find(yvec > errors.Y(1, 1) & yvec < errors.Y(end, 1));
        
    elseif dataset_index == 5
        % index limits for processing
        % X 320 to 392 (72); Y 205 to 325 (120)
        cmin_p = 5;
        cmax_p = 23; %303+18-1=320  303+90-1=392  X direction
        rmin_p = 5;
        rmax_p = 35; %188+18-1=205  188+138-1=325 Y direction
        
%         % indices to be used
%         r = rmin_p:rmax_p;
%         c = cmin_p:cmax_p;
        
        % only retain true solution in the index limits for
        % processing
        c = find(errors.X(1, :) >= 320  & errors.X(1, :) <= 392);
        r = find(errors.Y(:, 1) >= 205  & errors.Y(:, 1) <= 325);
    end

    grid_point_row_index = randsample(r, 1);
    grid_point_col_index = randsample(c, 1);

    %% extract errors for this grid point
    
    % identify indices on the true solution grid corresponding to the
    % current grid point
    if dataset_index <= 3
        r_t = grid_point_row_index;
        c_t = grid_point_col_index;
    elseif dataset_index == 4
        r_t = find(errors.Y(:, 1) == results.Y(grid_point_row_index, 1));
        c_t = find(errors.X(1, :) == results.X(1, grid_point_col_index));
    elseif dataset_index == 5
        % indices to be used in the true solution
%         rt = 8:2:69;
%         ct = 7:2:43;
%         r_t = grid_point_row_index - rmin_p + 1; %find(errors.Y(:, 1) == results.Y(grid_point_row_index, 1) + 187);
%         c_t = grid_point_col_index - cmin_p + 1; %find(errors.X(1, :) == results.X(1, grid_point_col_index) + 302);
        r_t = grid_point_row_index;
        c_t = grid_point_col_index;
    end
    
    % extract errors
    err_U_current = errors.err_U(r_t, c_t, snapshot_index);
    err_V_current = errors.err_V(r_t, c_t, snapshot_index);

    %% extract individual uncertainties for this grid point
    unc_current = struct;
    unc_current.imx = results.uncertainty2D.Uimx(grid_point_row_index, grid_point_col_index);
    unc_current.imy = results.uncertainty2D.Uimy(grid_point_row_index, grid_point_col_index);
    
    unc_current.mcx = results.uncertainty2D.MCx(grid_point_row_index, grid_point_col_index);
    unc_current.mcy = results.uncertainty2D.MCy(grid_point_row_index, grid_point_col_index);

    unc_current.csx = results.uncertainty2D.Ucsx(grid_point_row_index, grid_point_col_index);
    unc_current.csy = results.uncertainty2D.Ucsy(grid_point_row_index, grid_point_col_index);

    %% calculate local weights for this grid point
    
    % --------------------
    % load deformed image
    % --------------------
    im1 = imread(fullfile(files_im1{window_size_index, dataset_index}(snapshot_index).folder, files_im1{window_size_index, dataset_index}(snapshot_index).name));
    im1d = double(flipud(im1));
    
    im2 = imread(fullfile(files_im2{window_size_index, dataset_index}(snapshot_index).folder, files_im2{window_size_index, dataset_index}(snapshot_index).name));
    im2d = double(flipud(im2));
    
    % extract image height and width
    [image_height, image_width] = size(im1d);
    
    % -----------------------------
    % extract interrogation window
    % -----------------------------
    % extract grid-point image co-ordinates
    X = results.X(grid_point_row_index, grid_point_col_index);
    Y = results.Y(grid_point_row_index, grid_point_col_index);
    
    % extract window size
    str = strsplit(jobfile.Data.PIV4.winsize, ',');
    window_size_x = str2double(str{1});
    window_size_y = str2double(str{2});
    
    % extract window resolution
    str = strsplit(jobfile.Data.PIV4.winres, ',');
    window_resolution_x = str2double(str{1});
    window_resolution_y = str2double(str{1});
    
    % identify window extents
    xmin = X - window_size_x/2 + 1;
    xmax = X + window_size_x/2;
    ymin = Y - window_size_y/2 + 1;
    ymax = Y + window_size_y/2;
    
    %% extract image windows
    
    % find the image windows
    zone1 = im1d(max([1 ymin]):min([image_height ymax]), max([1 xmin]):min([image_width xmax]));
    zone2 = im2d(max([1 ymin]):min([image_height ymax]), max([1 xmin]):min([image_width xmax]));
    
    % minimum subtraction
    zone1 = zone1 - min(zone1(:));
    zone2 = zone2 - min(zone2(:));
    
    % zero pad the image windows
    if size(zone1,1)~=window_size_y || size(zone1,2)~=window_size_x
        w1 = zeros(window_size_y,window_size_x);
        w1( 1+max([0 1-ymin]):window_size_y - max([0 ymax-image_height]), 1+max([0 1-xmin]):window_size_x - max([0 xmax-image_width]) ) = zone1;
        zone1 = w1;
    end
    
    if size(zone2,1)~=window_size_y || size(zone2,2)~=window_size_x
        w2 = zeros(window_size_y,window_size_x);
        w2( 1+max([0 1-ymin]):window_size_y - max([0 ymax-image_height]), 1+max([0 1-xmin]):window_size_x -max([0 xmax-image_width]) ) = zone2;
        zone2 = w2;
    end
    
    if max(zone1(:)) == 0 || max(zone2(:)) == 0
        continue;
    end
    
    %% apply windowing function to interrogation window
    
    % create window masking filter
    sfilt1 = windowmask([window_size_x window_size_y], [window_resolution_x window_resolution_y]);
    sfilt2 = windowmask([window_size_x window_size_y], [window_resolution_x window_resolution_y]);
    
    % apply the image spatial filter
    im1_sub = zone1 .* sfilt1;
    im2_sub = zone2 .* sfilt2;
    
    %% identify particles
    
    % calculate product of images
    im_p = sqrt(im1_sub .* im2_sub); %/sqrt(max(im1_sub(:)) * max(im2_sub(:)));
    
    % set minimum intensity threshold for the image    
    if dataset_index == 3
        intensity_threshold = prctile(im_p(:), 99);
    elseif dataset_index == 5
        intensity_threshold = prctile(im_p(:), 95);
    else
        intensity_threshold = prctile(im_p(:), 97);
    end

    % identify intensity peaks using particle identification
    [p_matrix, peaks, num_p] = dynamic_threshold_segmentation_v3(im_p, intensity_threshold, 0);

    %% dot sizing and centroid estimation
    
    % extract cross-correlation diameter
    d_c = results.uncertainty2D.Autod(grid_point_row_index, grid_point_col_index);
    % calculate particle diameter from the auto-correlation
    d_p = d_c/sqrt(2);
    % change size properties
    % particle sizing settings
    sizeprops_current = sizeprops;
    sizeprops_current.p_area = 0.5*d_p;
    % perform particle sizing
    [XYDiameter, mapsizeinfo, locxy]=particle_size_MAIN_V1(im_p, p_matrix, num_p, sizeprops_current);
    
    % extract locations of particles
    x1 = XYDiameter(:, 1);
    y1 = XYDiameter(:, 2);
    
    % extract particle diameters
    d1 = XYDiameter(:, 3);
    
    % identify nan elements
    indices = isnan(x1) | isnan(y1);
    
    % remove nan elements
    x1(indices) = [];
    y1(indices) = [];
    d1(indices) = [];
    
    % assign the positions of the peaks in the second image
    % to be the same (as the disparity is within subpixel
    % limit)
    x2 = x1;
    y2 = y1;
    
    %% display identification results    
%     figure
%     imagesc(im_p)
%     colormap(gray)
%     caxis([0 intensity_threshold])
% %     maxval = 0.8 * max(im_p(:))/2;
% %     caxis([0 maxval])
%     hold on
%     plot(x1, y1, 'o')
%     
%     % plot window edge
%     x_c = size(im_p, 2)/2;
%     y_c = size(im_p, 1)/2;
%     
%     x = [x_c - window_resolution_x/2, x_c + window_resolution_x/2, x_c + window_resolution_x/2, x_c - window_resolution_x/2, x_c - window_resolution_x/2];
%     y = [y_c - window_resolution_y/2, y_c - window_resolution_y/2, y_c + window_resolution_y/2, y_c + window_resolution_y/2, y_c - window_resolution_y/2];
%     plot(x, y, '*-')
%     colorbar
%     set_axes(gca);
%     title(dataset_name_array{dataset_index});
%     pause(0.1);
    
    %% remove particles and perform resampling calc
    
    % calculate number of particles identified
    num_particles = numel(x1);
    % calculate number of particles to be removed
    num_particles_remove = round(percentage_particles_remove * num_particles);
    % calculate local variance
    % [unc_resampling_avg{trial_index}, unc_resampling_std{trial_index}, unc_resampling{trial_index}] = calculate_variance_resampling(im1_sub, im2_sub, window_size_x/2, window_size_y/2, 0, 0, x1, y1, x2, y2, d_p*ones(size(x1)), num_particles_remove, num_resampling_trials, uncertainty_flags, 0.25*intensity_threshold, [window_size_x, window_size_y], [window_resolution_x, window_resolution_y]);        
    [unc_resampling_avg{trial_index}, unc_resampling_std{trial_index}, unc_resampling{trial_index}] = calculate_variance_resampling(im1_sub, im2_sub, x1, y1, x2, y2, ...
    d_p*ones(size(x1)), num_particles_remove, num_resampling_trials, uncertainty_flags, 0.25*intensity_threshold, ...
    [window_size_x, window_size_y], [window_resolution_x, window_resolution_y]);

    %% store error and uncertainty
    
    unc_trials{trial_index} = unc_current;
    err_U(trial_index) = err_U_current;
    err_V(trial_index) = err_V_current;

end

%% save calculations to file

fprintf('saving results to file\n');
filename = fullfile(current_write_directory, 'monte_carlo_results.mat');
save(filename, 'err_U', 'err_V', 'unc_trials', ...
    'unc_resampling_avg', 'unc_resampling_std', 'unc_resampling');

%% plot histogram of datasets that have been accessed

dataset_count = nans(num_datasets, num_window_size);
X = categorical(dataset_name_array);
% find number of trials falling in each dataset
for dataset_index = 1:num_datasets
    for window_size_index = 1:num_window_size
        indices = find(dataset_index_array == dataset_index);
        dataset_count(dataset_index, window_size_index) = sum(window_resolution_index_array(indices) == window_size_index);
    end
end

figure
bar(X, dataset_count, 'stacked');
ylabel('Count')
legend('WR=32', 'WR=64', 'location', 'northoutside', 'Orientation', 'horizontal')

if save_figures
    save_figure_to_png_eps_fig(figure_write_directory, 'datasets-accessed', [true, false, false]);
end

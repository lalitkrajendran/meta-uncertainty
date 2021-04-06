clear
close all
clc

restoredefaultpath;
addpath(genpath('/scratch/shannon/c/aether/Projects/BOS/general-codes/matlab-codes/'));
addpath('../prana-uncertainty-average-dc-new-im-2-save-deform-cs/');
setup_default_settings;
% dbstop if error

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

% array of uncertainty combination model to be used
% 'unweighted', 'global-weight-var', 'global-weight-std',
% 'local-weight-var', 'local-weight-std'
uncertainty_combination_method_array = {'unweighted'; 'global-weight-var'; 'local-weight-var'};
num_uncertainty_combination_methods = numel(uncertainty_combination_method_array);

% directory where results of this analysis are to be saved
top_write_directory = fullfile('/scratch/shannon/a/aether/Projects/PIV_PTV_Uncertainty_Quantification/uncertainty-combination/analysis/results/', 'monte-carlo', 'new-processing-all-datasets');
mkdir_c(top_write_directory);

%% analysis settings

% number of trials
num_trials = 1e3;
% minimum allowable error (pix.)
min_error_threshold = 1e-4;
% maximum allowable error (pix.)
max_error_threshold = 0.2;
% number of bins for histogram
num_bins = 25;
% bins for histograms
bins = linspace(min_error_threshold, max_error_threshold, num_bins);
% number of bins for the coarse histogram
num_bins_coarse = 8;
% save figures? (True/False)
save_figures = true;

%% bootstrapping settings

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
sizeprops.p_area = 3; %0.5 * d_p^2;
sizeprops.sigma = 4;
sizeprops.errors = 0;
sizeprops.save_dir = [];

%% directory settings for this case

% directory to save results for this case
current_write_directory = fullfile(top_write_directory, ['max_error=' num2str(max_error_threshold, '%.2f') 'pix_trials=' num2str(num_trials, '%d')], ['ppr' num2str(percentage_particles_remove, '%.2f') '_nrs' num2str(num_resampling_trials)]);
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

%% calculate global uncertainty weights for all datasets

fprintf('Calculating global weights\n');

% loop through window resolutions
for window_size_index = 1:num_window_size
    % loop through datasets
    for dataset_index = 1:numel(dataset_name_array)
        % name of the data set
        dataset_name = dataset_name_array{dataset_index};
        fprintf('Dataset: %s\n', dataset_name);
        
        for snapshot_index = 1:numel(results_all{window_size_index, dataset_index})
            %% calculate weights

            results = results_all{window_size_index, dataset_index}(snapshot_index);
            
            %% select vector locations for which true solution is available
            if snapshot_index == 1
                if dataset_index <= 3
                   c = 1:size(results.X, 2);
                   r = 1:size(results.X, 1);
                elseif dataset_index == 4                
                    %% find vector locations in the processing results for which
                    % true solution is available

                    % extract co-ordinate arrays for the processing grid
                    xvec = results.X(1, :);
                    yvec = results.Y(:, 1);

                    % find vector locations in the processing results for which
                    % true solution is available
                    c = find(xvec > errors_all{window_size_index, dataset_index}.X(1, 1) & xvec < errors_all{window_size_index, dataset_index}.X(1, end));
                    r = find(yvec > errors_all{window_size_index, dataset_index}.Y(1, 1) & yvec < errors_all{window_size_index, dataset_index}.Y(end, 1));
                elseif dataset_index == 5
                    % index limits for processing
                    % X 320 to 392 (72); Y 205 to 325 (120)
                    cmin_p = 5;
                    cmax_p = 23; %303+18-1=320  303+90-1=392  X direction
                    rmin_p = 5;
                    rmax_p = 35; %188+18-1=205  188+138-1=325 Y direction
                    
                    % indices to be used
                    r = rmin_p:rmax_p;
                    c = cmin_p:cmax_p;
                end
            end
            
            %% calculate weights
            % assign weights as the inverse of the variance of the
            % uncertainty in field of view

            % IM
            results_all{window_size_index, dataset_index}(snapshot_index).w_glob_im_x = 1/var(results.uncertainty2D.Uimx(r, c), [], 'all', 'omitnan');
            results_all{window_size_index, dataset_index}(snapshot_index).w_glob_im_y = 1/var(results.uncertainty2D.Uimy(r, c), [], 'all', 'omitnan');
            
            % MC
            results_all{window_size_index, dataset_index}(snapshot_index).w_glob_mc_x = 1/var(results.uncertainty2D.MCx(r, c), [], 'all', 'omitnan');
            results_all{window_size_index, dataset_index}(snapshot_index).w_glob_mc_y = 1/var(results.uncertainty2D.MCy(r, c), [], 'all', 'omitnan');
            
            % CS
            results_all{window_size_index, dataset_index}(snapshot_index).w_glob_cs_x = 1/var(real(results.uncertainty2D.Ucsx(r, c)), [], 'all', 'omitnan');
            results_all{window_size_index, dataset_index}(snapshot_index).w_glob_cs_y = 1/var(real(results.uncertainty2D.Ucsy(r, c)), [], 'all', 'omitnan');
        end 
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
sigma_im_x = nans(1, num_trials);
sigma_im_y = nans(1, num_trials);
sigma_mc_x = nans(1, num_trials);
sigma_mc_y = nans(1, num_trials);
sigma_cs_x = nans(1, num_trials);
sigma_cs_y = nans(1, num_trials);

sigma_unwt_x = nans(1, num_trials);
sigma_unwt_y= nans(1, num_trials);
sigma_glob_x = nans(1, num_trials);
sigma_glob_y= nans(1, num_trials);
sigma_loc_x = nans(1, num_trials);
sigma_loc_y= nans(1, num_trials);

% global weights
w_glob_im_x = nans(1, num_trials);
w_glob_im_y = nans(1, num_trials);
w_glob_mc_x = nans(1, num_trials);
w_glob_mc_y = nans(1, num_trials);
w_glob_cs_x = nans(1, num_trials);
w_glob_cs_y = nans(1, num_trials);

% local weights
w_loc_im_x = nans(1, num_trials);
w_loc_im_y = nans(1, num_trials);
w_loc_mc_x = nans(1, num_trials);
w_loc_mc_y = nans(1, num_trials);
w_loc_cs_x = nans(1, num_trials);
w_loc_cs_y = nans(1, num_trials);

% uncertainty from resampling
sigma_trials = cell(1, num_trials);
%% loop through trials and perform analysis

fprintf('running monte-carlo\n');

% set seed for random number generator
rng(0);

% generate random numbers for overall dataset index
% overall_dataset_index_array = 1 + round(rand(1, num_trials) * (num_datasets * num_window_resolution - 1));
overall_dataset_index_array = 1 + round(rand(1, num_trials) * (num_datasets - 1));

overall_window_resolution_index_array = 1 + round(rand(1, num_trials) * (num_window_size - 1));

for trial_index = 1:num_trials
    fprintf('trial_index: %d\n', trial_index);

    % calculate random number for data set
    dataset_index_array(trial_index) = randi(num_datasets);
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
       
    %% extract weights for this snapshot
        
    % assign weights as the inverse of the variance of the
    % uncertainty in field of view
    
    % IM
    w_glob.im_x(trial_index) = results.w_glob_im_x;
    w_glob_im_y(trial_index) = results.w_glob_im_y;
    
    % MC
    w_glob_mc_x(trial_index) = results.w_glob_mc_x;
    w_glob_mc_y(trial_index) = results.w_glob_mc_y;

    % CS
    w_glob_cs_x(trial_index) = results.w_glob_cs_x;
    w_glob_cs_y(trial_index) = results.w_glob_cs_y;

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
        
        % indices to be used
        r = rmin_p:rmax_p;
        c = cmin_p:cmax_p;
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
        rt = 8:2:69;
        ct = 7:2:43;
        r_t = grid_point_row_index - rmin_p + 1; %find(errors.Y(:, 1) == results.Y(grid_point_row_index, 1) + 187);
        c_t = grid_point_col_index - cmin_p + 1; %find(errors.X(1, :) == results.X(1, grid_point_col_index) + 302);
    end
    
    % extract errors
    err_U(trial_index) = errors.err_U(r_t, c_t, snapshot_index);
    err_V(trial_index) = errors.err_V(r_t, c_t, snapshot_index);

    %% extract individual uncertainties for this grid point
    
    sigma_im_x(trial_index) = results.uncertainty2D.Uimx(grid_point_row_index, grid_point_col_index);
    sigma_im_y(trial_index) = results.uncertainty2D.Uimy(grid_point_row_index, grid_point_col_index);

    sigma_mc_x(trial_index) = results.uncertainty2D.MCx(grid_point_row_index, grid_point_col_index);
    sigma_mc_y(trial_index) = results.uncertainty2D.MCy(grid_point_row_index, grid_point_col_index);

    sigma_cs_x(trial_index) = real(results.uncertainty2D.Ucsx(grid_point_row_index, grid_point_col_index));
    sigma_cs_y(trial_index) = real(results.uncertainty2D.Ucsy(grid_point_row_index, grid_point_col_index));

    %% calculate combined uncertainty for that point for unweighted average
    
    sigma_unwt_x(trial_index) = 1/3 .* (sigma_im_x(trial_index) + sigma_mc_x(trial_index) + sigma_cs_x(trial_index));
    sigma_unwt_y(trial_index) = 1/3 .* (sigma_im_y(trial_index) + sigma_mc_y(trial_index) + sigma_cs_y(trial_index));
    
    %% calculate combined uncertainty for that point for global weighted average
    
    sigma_glob_x(trial_index) = 1./(w_glob.im_x(trial_index) + w_glob_mc_x(trial_index) + w_glob_cs_x(trial_index)) .* (w_glob.im_x(trial_index) .* sigma_im_x(trial_index) + w_glob_mc_x(trial_index) .* sigma_mc_x(trial_index) + w_glob_cs_x(trial_index) .* sigma_cs_x(trial_index));
    sigma_glob_y(trial_index) = 1./(w_glob_im_y(trial_index) + w_glob_mc_y(trial_index) + w_glob_cs_y(trial_index)) .* (w_glob_im_y(trial_index) .* sigma_im_y(trial_index) + w_glob_mc_y(trial_index) .* sigma_mc_y(trial_index) + w_glob_cs_y(trial_index) .* sigma_cs_y(trial_index));

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
    zone1 = im1d( max([1 ymin]):min([image_height ymax]),max([1 xmin]):min([image_width xmax]));
    zone2 = im2d( max([1 ymin]):min([image_height ymax]),max([1 xmin]):min([image_width xmax]));
    
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
    
    % calculate minimum intensity threshold for the image
    intensity_threshold = min(im1d(:));
    
    % calculate product of images
    im_p = im1_sub.*im2_sub/max(im1_sub(:));
    
    % identify intensity peaks using particle identification
    [p_matrix,peaks,num_p] = dynamic_threshold_segmentation_v3(im_p, intensity_threshold, 0);
    
    % perform particle sizing
    [XYDiameter, mapsizeinfo, locxy]=particle_size_MAIN_V1(im_p, p_matrix, num_p, sizeprops);
    
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
    
    %% remove particles and perform resampling calc
    
    % calculate number of particles identified
    num_particles = numel(x1);
    % calculate number of particles to be removed
    num_particles_remove = round(percentage_particles_remove * num_particles);
    % calculate local variance
    [sigma_avg, sigma_std, sigma_trials{trial_index}] = calculate_variance_resampling(im1_sub, im2_sub, window_size_x/2, window_size_y/2, 0, 0, x1, y1, x2, y2, d_p*ones(size(x1)), num_particles_remove, num_resampling_trials, uncertainty_flags, intensity_threshold, [window_size_x, window_size_y], [window_resolution_x, window_resolution_y]);
    
%     % temporary fix
%     if sigma_std.mcx > 5 * sigma_std.mcy
%         sigma_std.mcx = sigma_std.mcy;
%     end
%     if sigma_std.csy > 5 * sigma_std.csx
%         sigma_std.csy = sigma_std.csx;
%     end
    
    %% calculate local weights

    w_loc_im_x(trial_index) = 1/sigma_std.imx^2/(1/sigma_std.imx^2 + 1/sigma_std.mcx^2 + 1/sigma_std.csx^2);
    w_loc_im_y(trial_index) = 1/sigma_std.imy^2/(1/sigma_std.imy^2 + 1/sigma_std.mcy^2 + 1/sigma_std.csy^2);

    w_loc_mc_x(trial_index) = 1/sigma_std.mcx^2/(1/sigma_std.imx^2 + 1/sigma_std.mcx^2 + 1/sigma_std.csx^2);
    w_loc_mc_y(trial_index) = 1/sigma_std.mcy^2/(1/sigma_std.imy^2 + 1/sigma_std.mcy^2 + 1/sigma_std.csy^2);

    w_loc_cs_x(trial_index) = 1 - (w_loc_im_x(trial_index) + w_loc_mc_x(trial_index));
    w_loc_cs_y(trial_index) = 1 - (w_loc_im_y(trial_index) + w_loc_mc_y(trial_index));

    %% calculate locally weighted uncertainty    
    
    sigma_loc_x(trial_index) = 1./(w_loc_im_x(trial_index) + w_loc_mc_x(trial_index) + w_loc_cs_x(trial_index)) .* (w_loc_im_x(trial_index) .* sigma_im_x(trial_index) + w_loc_mc_x(trial_index) .* sigma_mc_x(trial_index) + w_loc_cs_x(trial_index) .* sigma_cs_x(trial_index));
    sigma_loc_y(trial_index) = 1./(w_loc_im_y(trial_index) + w_loc_mc_y(trial_index) + w_loc_cs_y(trial_index)) .* (w_loc_im_y(trial_index) .* sigma_im_y(trial_index) + w_loc_mc_y(trial_index) .* sigma_mc_y(trial_index) + w_loc_cs_y(trial_index) .* sigma_cs_y(trial_index));
    
    
end

%% save calculations to file
save(filename, '
%% filter out invalid measurements

% remove invalid measurements for the error
err_all_valid = nan_invalid_measurements([err_U(:); err_V(:)], min_error_threshold, max_error_threshold);

% remove invalid measurements for IM uncertainty
sigma_im_all_valid = nan_invalid_measurements([sigma_im_x(:); sigma_im_y(:)], min_error_threshold, max_error_threshold);
% remove invalid measurements for MC uncertainty
sigma_mc_all_valid = nan_invalid_measurements([sigma_mc_x(:); sigma_mc_y(:)], min_error_threshold, max_error_threshold);
% remove invalid measurements for CS uncertainty
sigma_cs_all_valid = nan_invalid_measurements([sigma_cs_x(:); sigma_cs_y(:)], min_error_threshold, max_error_threshold);
% remove invalid measurements for the unweighted uncertainty
sigma_unwt_all_valid = nan_invalid_measurements([sigma_unwt_x(:); sigma_unwt_y(:)], min_error_threshold, max_error_threshold);
% remove invalid measurements for the globally weighted uncertainty
sigma_glob_all_valid = nan_invalid_measurements([sigma_glob_x(:); sigma_glob_y(:)], min_error_threshold, max_error_threshold);
% remove invalid measurements for the locally weighted uncertainty
sigma_loc_all_valid = nan_invalid_measurements([sigma_loc_x(:); sigma_loc_y(:)], min_error_threshold, max_error_threshold);

%% calculate statistics

% calculate rms of error
err_rms = rms(err_all_valid, 'omitnan');

% calculate rms of IM uncertainty
sigma_im_rms = rms(sigma_im_all_valid, 'omitnan');
% calculate rms of MC uncertainty
sigma_mc_rms = rms(sigma_mc_all_valid, 'omitnan');
% calculate rms of CS uncertainty
sigma_cs_rms = rms(sigma_cs_all_valid, 'omitnan');
% calculate rms of unweighted uncertainty
sigma_unwt_rms = rms(sigma_unwt_all_valid, 'omitnan');
% calculate rms of globally weighted uncertainty
sigma_glob_rms = rms(sigma_glob_all_valid, 'omitnan');
% calculate rms of locally weighted uncertainty
sigma_loc_rms = rms(sigma_loc_all_valid, 'omitnan');

%% calculate pdf of the error and uncertainty distributions

% calculate pdf of error
[N_err, ~] = histcounts(err_all_valid, bins, 'normalization', 'pdf');

% calculate pdf of IM uncertainty
[N_sigma_im, ~] = histcounts(sigma_im_all_valid, bins, 'normalization', 'pdf');
% calculate pdf of MC uncertainty
[N_sigma_mc, ~] = histcounts(sigma_mc_all_valid, bins, 'normalization', 'pdf');
% calculate pdf of CS uncertainty
[N_sigma_cs, ~] = histcounts(sigma_cs_all_valid, bins, 'normalization', 'pdf');
% calculate pdf of unweighted uncertainty
[N_sigma_unwt, ~] = histcounts(sigma_unwt_all_valid, bins, 'normalization', 'pdf');
% calculate pdf of globally weighted uncertainty
[N_sigma_glob, ~] = histcounts(sigma_glob_all_valid, bins, 'normalization', 'pdf');
% calculate pdf of locally weighted uncertainty
[N_sigma_loc, ~] = histcounts(sigma_loc_all_valid, bins, 'normalization', 'pdf');

%% calculate coverage

% true coverage
cov_true = calculate_coverage_percentage(err_all_valid, err_rms);
% im coverage
cov_im = calculate_coverage_percentage(err_all_valid, sigma_im_all_valid);
% mc coverage
cov_mc = calculate_coverage_percentage(err_all_valid, sigma_mc_all_valid);
% cs coverage
cov_cs = calculate_coverage_percentage(err_all_valid, sigma_cs_all_valid);
% unweighted coverage
cov_unwt = calculate_coverage_percentage(err_all_valid, sigma_unwt_all_valid);
% globally weighted coverage
cov_glob = calculate_coverage_percentage(err_all_valid, sigma_glob_all_valid);
% locally weighted coverage
cov_loc = calculate_coverage_percentage(err_all_valid, sigma_loc_all_valid);

%% calculate rms error and uncertainty binwise

% calculate binwise rms of error and uncertainty for im
[err_rms_binwise_im, sigma_rms_binwise_im] = calculate_rms_binwise(err_all_valid, sigma_im_all_valid, max_error_threshold, num_bins_coarse);
% calculate binwise rms of error and uncertainty for mc
[err_rms_binwise_mc, sigma_rms_binwise_mc] = calculate_rms_binwise(err_all_valid, sigma_mc_all_valid, max_error_threshold, num_bins_coarse);
% calculate binwise rms of error and uncertainty for cs
[err_rms_binwise_cs, sigma_rms_binwise_cs] = calculate_rms_binwise(err_all_valid, sigma_cs_all_valid, max_error_threshold, num_bins_coarse);
% calculate binwise rms of error and uncertainty for unweighted
[err_rms_binwise_unwt, sigma_rms_binwise_unwt] = calculate_rms_binwise(err_all_valid, sigma_unwt_all_valid, max_error_threshold, num_bins_coarse);
% calculate binwise rms of error and uncertainty for globally weighted
[err_rms_binwise_glob, sigma_rms_binwise_glob] = calculate_rms_binwise(err_all_valid, sigma_glob_all_valid, max_error_threshold, num_bins_coarse);
% calculate binwise rms of error and uncertainty for locally weighted
[err_rms_binwise_loc, sigma_rms_binwise_loc] = calculate_rms_binwise(err_all_valid, sigma_loc_all_valid, max_error_threshold, num_bins_coarse);

%% calculate histogram of global weights

% bins for histogram
bins_weights = linspace(0, 1, 25);

% calculate normalized weights for all methods, combining x and y weights
w_im_all_glob = [w_glob.im_x./(w_glob.im_x + w_glob_mc_x + w_glob_cs_x), ...
    w_glob_im_y./(w_glob_im_y + w_glob_mc_y + w_glob_cs_y)];
w_mc_all_glob = [w_glob_mc_x./(w_glob.im_x + w_glob_mc_x + w_glob_cs_x), ...
    w_glob_mc_y./(w_glob_im_y + w_glob_mc_y + w_glob_cs_y)];
w_cs_all_glob = [w_glob_cs_x./(w_glob.im_x + w_glob_mc_x + w_glob_cs_x), ...
    w_glob_cs_y./(w_glob_im_y + w_glob_mc_y + w_glob_cs_y)];

% calculate pdf of image matching weights
[N_w_im_glob, ~] = histcounts(w_im_all_glob, bins_weights, 'normalization', 'pdf');
% calculate pdf of mc weights
[N_w_mc_glob, ~] = histcounts(w_mc_all_glob, bins_weights, 'normalization', 'pdf');
% calculate pdf of cs weights
[N_w_cs_glob, ~] = histcounts(w_cs_all_glob, bins_weights, 'normalization', 'pdf');

%% calculate histogram of local weights

% bins for histogram
bins_weights = linspace(0, 1, 25);

% calculate normalized weights for all methods, combining x and y weights
w_im_all_loc = [w_loc_im_x./(w_loc_im_x + w_loc_mc_x + w_loc_cs_x), ...
    w_loc_im_y./(w_loc_im_y + w_loc_mc_y + w_loc_cs_y)];
w_mc_all_loc = [w_loc_mc_x./(w_loc_im_x + w_loc_mc_x + w_loc_cs_x), ...
    w_loc_mc_y./(w_loc_im_y + w_loc_mc_y + w_loc_cs_y)];
w_cs_all_loc = [w_loc_cs_x./(w_loc_im_x + w_loc_mc_x + w_loc_cs_x), ...
    w_loc_cs_y./(w_loc_im_y + w_loc_mc_y + w_loc_cs_y)];

% calculate pdf of image matching weights
[N_w_im_loc, ~] = histcounts(w_im_all_loc, bins_weights, 'normalization', 'pdf');
% calculate pdf of mc weights
[N_w_mc_loc, ~] = histcounts(w_mc_all_loc, bins_weights, 'normalization', 'pdf');
% calculate pdf of cs weights
[N_w_cs_loc, ~] = histcounts(w_cs_all_loc, bins_weights, 'normalization', 'pdf');

%% write rms values to file

fprintf('writing statistics to file\n');

fileID = fopen(fullfile(current_write_directory, 'rms-values.txt'), 'w');

fprintf(fileID, '-------------------------------\n');
fprintf(fileID, 'Errors (pix.)\n');
fprintf(fileID, '-------------------------------\n');
fprintf(fileID, 'Errors: %.3f\n', err_rms);

fprintf(fileID, '-------------------------------\n');
fprintf(fileID, 'Uncertainties (pix.)\n');
fprintf(fileID, '-------------------------------\n');

fprintf(fileID, 'IM: %.3f\n', sigma_im_rms);
fprintf(fileID, 'MC: %.3f\n', sigma_mc_rms);
fprintf(fileID, 'CS: %.3f\n', sigma_cs_rms);
fprintf(fileID, 'Unweighted: %.3f\n', sigma_unwt_rms);
fprintf(fileID, 'Global: %.3f\n', sigma_glob_rms);
fprintf(fileID, 'Local: %.3f\n', sigma_loc_rms);

fclose(fileID);

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

%% plot histogram of weights

% ====================================
% histogram of normalized weights
% ====================================
figure

bins_plot = bins_weights(1:end-1);
y_max = 10;

% IM
l_im_glob = plot(bins_plot, N_w_im_glob, 'c');
hold on
l_im_loc = plot(bins_plot, N_w_im_loc, 'c--');
% MC
l_mc_glob = plot(bins_plot, N_w_mc_glob, 'r');
l_mc_loc = plot(bins_plot, N_w_mc_loc, 'r--');
% CS
l_cs_glob = plot(bins_plot, N_w_cs_glob, 'm');
l_cs_loc = plot(bins_plot, N_w_cs_loc, 'm--');

ylim([0 y_max])
xlabel('Weights')
ylabel('Probability Density')
legend([l_im_glob, l_im_loc, l_mc_glob, l_mc_loc, l_cs_glob, l_cs_loc], ...
    {'IM, global', 'IM, local', 'MC, global', 'MC, local', 'CS, global', 'CS, local'}, ...
    'location', 'eastoutside')
set(gcf, 'Position', [336   482   583   412])

if save_figures
    save_figure_to_png_eps_fig(figure_write_directory, 'normalized-weights-histogram', [true, false, false]);
end

%% plot histograms of errors and uncertainties with rms 

% =================================
% Error and Uncertainty histograms
% =================================

bins_plot = bins(1:end-1);
y_max = 30;

figure
% average error
l_err_comb = plot(bins_plot, N_err, 'k');
hold on
plot(err_rms*[1 1], [0 y_max], 'k');
% IM
l_im = plot(bins_plot, N_sigma_im, 'c');
plot(sigma_im_rms * [1, 1], [0 y_max], 'c')
% MC
l_mc = plot(bins_plot, N_sigma_mc, 'r');
plot(sigma_mc_rms * [1, 1], [0 y_max], 'r')
% CS
l_cs = plot(bins_plot, N_sigma_cs, 'm');
plot(sigma_cs_rms * [1, 1], [0 y_max], 'm')
% unweighted uncertainty
l_unwt = plot(bins_plot, N_sigma_unwt, 'g');
plot(sigma_unwt_rms*[1 1], [0 y_max], 'g')
% globally weighted uncertainty
l_glob = plot(bins_plot, N_sigma_glob, 'g--');
plot(sigma_glob_rms*[1 1], [0 y_max], 'g--')
% locally weighted uncertainty
l_loc = plot(bins_plot, N_sigma_loc, 'g:');
plot(sigma_loc_rms*[1 1], [0 y_max], 'g:')

ylim([0 y_max])
xlabel('(pix.)')
ylabel('Probability Density')
legend([l_err_comb, l_im, l_mc, l_cs, l_unwt, l_glob, l_loc], ...
{'Error', '\sigma_{IM}', '\sigma_{MC}', '\sigma_{CS}', '\sigma_{unwt}', '\sigma_{glob}', '\sigma_{loc}'}, ...
'location', 'eastoutside')
set(gcf, 'Position', [336   482   583   412])

if save_figures
    save_figure_to_png_eps_fig(figure_write_directory, 'histograms', [true, false, false]);
end

%% plot rms error vs rms uncertainty

% ==============================
% RMS Error vs RMS Uncertainty
% ==============================
figure
% Error
plot([0 max_error_threshold], [0 max_error_threshold], 'k--');
hold on
% IM
l_im = plot(err_rms_binwise_im, sigma_rms_binwise_im, 'co-', 'markerfacecolor', 'c');
% MC
l_mc = plot(err_rms_binwise_mc, sigma_rms_binwise_mc, 'ro-', 'markerfacecolor', 'r');
% CS
l_cs = plot(err_rms_binwise_cs, sigma_rms_binwise_cs, 'mo-', 'markerfacecolor', 'm');
% unweighted uncertainty
l_unwt = plot(err_rms_binwise_unwt, sigma_rms_binwise_unwt, 'go-');
% globally weighted uncertainty
l_glob = plot(err_rms_binwise_glob, sigma_rms_binwise_glob, 'g^-');
% locally weighted uncertainty
l_loc = plot(err_rms_binwise_loc, sigma_rms_binwise_loc, 'gs-');

grid on
axis equal
axis([0 max_error_threshold 0 max_error_threshold])
xlabel('RMS Error (pix.)')
ylabel('RMS Uncertainty (pix.)')
legend([l_im, l_mc, l_cs, l_unwt, l_glob, l_loc], ...
    {'IM', 'MC', 'CS', 'Unwt', 'Global', 'Local'}, 'location', 'eastoutside')
set(gcf, 'Position', [336   511   488   383])

if save_figures
    save_figure_to_png_eps_fig(figure_write_directory, 'rms-error-uncertainty', [true, false, false]);
end

%% plot coverage

% ====================================
% coverage
% ====================================

X = categorical({'IM', 'MC', 'CS', 'Unwt', 'Glob', 'Local'});
% sort in the desired order
X = reordercats(X,{'IM', 'MC', 'CS', 'Unwt', 'Glob', 'Local'});
Y = [cov_im, cov_mc, cov_cs, cov_unwt, cov_glob, cov_loc];

figure
bar(X,Y); %,'FaceColor','flat');
hold on
l_true = plot(X, cov_true * [1, 1, 1, 1, 1, 1], 'r', 'linewidth', 3.0);

ylim([0 100])
ylabel('Coverage (%)')
% legend(l_true, 'Target', 'location', 'northoutside')
set(gcf, 'Position', [336   482   583   412])
if save_figures
    save_figure_to_png_eps_fig(figure_write_directory, 'coverage', [true, false, false]);
end

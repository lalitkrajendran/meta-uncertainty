clear
close all
clc

addpath(genpath('/scratch/shannon/c/aether/Projects/BOS/general-codes/matlab-codes/'));
addpath('prana-uncertainty-average-dc-new-im-2/');
setup_default_settings;
dbstop if error
%% read/write settings

% window resolution
window_resolution = 64;
% window size
window_size = 2 * window_resolution;

% directory containing images
top_image_directory = fullfile('/scratch/shannon/a/aether/Projects/PIV_PTV_Uncertainty_Quantification/Planar_uncertainty_work/Images/');
% directory containing files to be read
top_read_dirsectory = fullfile('/scratch/shannon/a/aether/Projects/PIV_PTV_Uncertainty_Quantification/Planar_uncertainty_work/Results/07_30_2018/', ['WS' num2str(window_resolution)], '/SOC_gradient/matfiles/');

% array of data set names
dataset_name_array = {'PivChal03B'; 'PivChal05B'; 'stagnation_flow'; 'Vortex_Ring'; 'Jetdata'};

% array of uncertainty combination model to be used
% 'unweighted', 'global-weight-var', 'global-weight-std',
% 'local-weight-var', 'local-weight-std'
uncertainty_combination_method_array = {'unweighted'; 'global-weight-var'; 'local-weight-var'};
num_uncertainty_combination_methods = numel(uncertainty_combination_method_array);

% directory where results of this analysis are to be saved
top_write_directory = fullfile('/scratch/shannon/a/aether/Projects/PIV_PTV_Uncertainty_Quantification/uncertainty-combination/analysis/results/', ['WS' num2str(window_resolution)]);
if ~exist(top_write_directory, 'dir')
    mkdir(top_write_directory);
end

%% plot settings

% save figures? (True/False)
save_figures = true;

% minimum uncertainty for a measurement to be considered valid
min_error_threshold = 1e-4;

% range of displacements to be displayed in the contour plots
displacement_color_min = 0;
displacement_color_max = 15;
displacement_contour_levels = linspace(displacement_color_min, displacement_color_max, 100);

% number of bins for the coarse histogram
num_bins = 8;

% uncertainty flags for the various methods (0 = False, 1 = True)
uncertainty_flags.ppruncertainty = 0;
uncertainty_flags.miuncertainty = 0;
uncertainty_flags.mcuncertainty = 1;

%% bootstrapping settings

% method to obtain particle position ('true', 'id')
particle_position_calc_method = 'id';
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

% number of particles to remove
% num_particles_remove_array = [5, 10, 15];
percentage_particles_remove = 0.2;
% number of resampling trials
num_resampling_trials = 10;
%% perform analysis
% data set that is to be loaded
for dataset_index = [1] %1:numel(dataset_name_array)
    % name of the data set
    dataset_name = dataset_name_array{dataset_index};
    fprintf('Dataset: %s\n', dataset_name);
    
    % directory containing images for this dataset
    current_image_directory = fullfile(top_image_directory, dataset_name);
    % get list of images in the directory
    [files, num_files] = get_directory_listing(current_image_directory, '*.tif', 'staticmask.tif');
    %% Load data

    % load results for vectors, errors and uncertainties
    results = load(fullfile(top_read_directory, [dataset_name '.mat']));
    % load results for histograms
    results_hist = load(fullfile(top_read_directory, [dataset_name '_hist.mat']));
    % load results for rms values of the errors and uncertainties
    results_rms = load(fullfile(top_read_directory, [dataset_name '_RMS.mat']));
    % load coverage results
    results_coverage = load(fullfile(top_read_directory, 'Coverage_07_30_2018_lowcutoff.mat'));

    % calculate size of the results array
    [num_rows, num_cols, num_snapshots] = size(results.Up);

    % account for peculiarities of each data set
    if dataset_index == 2
        results.Xd = fliplr(results.Xd');
        results.Yd = results.Yd';
    end    

    %% set uncertainty range for the data set
    if dataset_index == 3 || dataset_index == 5
        % maximum error for a measurement to be considered valid
        max_error_threshold = 0.3;

        % range of uncertainties to be displayed (pix.)
        uncertainty_color_min = 0;
        uncertainty_color_max = 0.3;
        uncertainty_contour_levels = linspace(uncertainty_color_min, uncertainty_color_max, 100);
        max_bin = 0.3;
    else
        % maximum error for a measurement to be considered valid
        max_error_threshold = 0.1;

        % range of uncertainties to be displayed (pix.)
        uncertainty_color_min = 0;
        uncertainty_color_max = 0.1;
        uncertainty_contour_levels = linspace(uncertainty_color_min, uncertainty_color_max, 100);
        max_bin = 0.08;
    end
    
    %%
    % ===================================================
    % loop through all uncertainty combination models
    % ===================================================
    
    for uncertainty_combination_method_index = 3 %1:numel(uncertainty_combination_method_array)
        % extract name of current uncertainty combination method
        uncertainty_combination_method = uncertainty_combination_method_array{uncertainty_combination_method_index};
        fprintf('method: %s\n', uncertainty_combination_method);
        
        % directory where current results are to be stored
        current_write_directory = fullfile(top_write_directory, dataset_name, uncertainty_combination_method);
        mkdir(current_write_directory);

        % directory to save figures
        figure_write_directory = fullfile(current_write_directory, 'figures');
        mkdir(figure_write_directory);

        %% allocate arrays
        U_avg = nans(num_rows, num_cols);
        V_avg = nans(num_rows, num_cols);

        sigma_U_combined = nans(num_rows, num_cols);
        sigma_V_combined = nans(num_rows, num_cols);

        err_U_combined = nans(num_rows, num_cols);
        err_V_combined = nans(num_rows, num_cols);

        w_im_x = nans(num_rows, num_cols, num_snapshots);
        w_im_y = nans(num_rows, num_cols, num_snapshots);

        w_mc_x = nans(num_rows, num_cols, num_snapshots);
        w_mc_y = nans(num_rows, num_cols, num_snapshots);

        w_cs_x = nans(num_rows, num_cols, num_snapshots);
        w_cs_y = nans(num_rows, num_cols, num_snapshots);

        %%
        % ===================================================
        % loop through snapshots and calculate average values
        % ===================================================

        for snapshot_index = 1:num_snapshots
            fprintf('snapshot_index: %d\n', snapshot_index);
            
            % interpolate davis displacement results onto prana grid
            Ud_interp = interp2(results.Xd, results.Yd, results.Ud(:, :, snapshot_index), results.Xp, results.Yp, 'cubic', 0);
            Vd_interp = interp2(results.Xd, results.Yd, results.Vd(:, :, snapshot_index), results.Xp, results.Yp, 'cubic', 0);

            % interpolate davis displacement uncertainty results onto prana grid
            sigma_Ud_interp = interp2(results.Xd, results.Yd, results.UCSx(:, :, snapshot_index), results.Xp, results.Yp, 'cubic', 0);
            sigma_Vd_interp = interp2(results.Xd, results.Yd, results.UCSy(:, :, snapshot_index), results.Xp, results.Yp, 'cubic', 0);

            % interpolate davis errors on to the prana grid
            err_Ud_interp = interp2(results.Xd, results.Yd, results.err_ud(:, :, snapshot_index), results.Xp, results.Yp, 'cubic', 0);
            err_Vd_interp = interp2(results.Xd, results.Yd, results.err_vd(:, :, snapshot_index), results.Xp, results.Yp, 'cubic', 0);

            % calculate average displacement
            U_avg(:, :, snapshot_index) = 1/2 * (results.Up(:, :, snapshot_index) + Ud_interp);
            V_avg(:, :, snapshot_index) = 1/2 * (results.Vp(:, :, snapshot_index) + Vd_interp);

            %% calculate weights for averaging
            if strcmpi(uncertainty_combination_method, 'unweighted')
                % assign equal weights for all methods

                % IM
                w_im_x(:, :, snapshot_index) = 1;
                w_im_y(:, :, snapshot_index) = 1;
                % MC
                w_mc_x(:, :, snapshot_index) = 1;
                w_mc_y(:, :, snapshot_index) = 1;
                % CS
                w_cs_x(:, :, snapshot_index) = 1;
                w_cs_y(:, :, snapshot_index) = 1;               
            elseif contains(uncertainty_combination_method, 'global-weight')
                if contains(uncertainty_combination_method, 'var')
                    % assign weights as the inverse of the variance of the
                    % uncertainty in field of view

                    % IM
                    w_im_x(:, :, snapshot_index) = 1/var(reshape(results.UIMx(:, :, snapshot_index), num_rows*num_cols, 1), [], 'omitnan');
                    w_im_y(:, :, snapshot_index) = 1/var(reshape(results.UIMy(:, :, snapshot_index), num_rows*num_cols, 1), [], 'omitnan');
                    % MC
                    w_mc_x(:, :, snapshot_index) = 1/var(reshape(results.UMCx(:, :, snapshot_index), num_rows*num_cols, 1), [], 'omitnan');
                    w_mc_y(:, :, snapshot_index) = 1/var(reshape(results.UMCy(:, :, snapshot_index), num_rows*num_cols, 1), [], 'omitnan');
                    % CS
                    w_cs_x(:, :, snapshot_index) = 1/var(sigma_Ud_interp(:), [], 'omitnan');
                    w_cs_y(:, :, snapshot_index) = 1/var(sigma_Vd_interp(:), [], 'omitnan');
                elseif contains(uncertainty_combination_method, 'std')
                    % assign weights as the inverse of the standard deviation of 
                    % the uncertainty in field of view

                    % IM
                    w_im_x(:, :, snapshot_index) = 1/std(reshape(results.UIMx(:, :, snapshot_index), num_rows*num_cols, 1), [], 'omitnan');
                    w_im_y(:, :, snapshot_index) = 1/std(reshape(results.UIMy(:, :, snapshot_index), num_rows*num_cols, 1), [], 'omitnan');
                    % MC
                    w_mc_x(:, :, snapshot_index) = 1/std(reshape(results.UMCx(:, :, snapshot_index), num_rows*num_cols, 1), [], 'omitnan');
                    w_mc_y(:, :, snapshot_index) = 1/std(reshape(results.UMCy(:, :, snapshot_index), num_rows*num_cols, 1), [], 'omitnan');
                    % CS
                    w_cs_x(:, :, snapshot_index) = 1/std(sigma_Ud_interp(:), [], 'omitnan');
                    w_cs_y(:, :, snapshot_index) = 1/std(sigma_Vd_interp(:), [], 'omitnan');
                end
            elseif contains(uncertainty_combination_method, 'local-weight')
                %% load image pair
                im1 = imread(fullfile(current_image_directory, files(2*(snapshot_index-1) + 1).name));
                im2 = imread(fullfile(current_image_directory, files(2*(snapshot_index)).name));
                
                % covert to double
                im1 = double(im1);
                im2 = double(im2);
                
                % calculate image size
                [image_height, image_width] = size(im1);
                
                %% extract displacement field for this snapshot
                X = results.Xp;
                Y = results.Yp;
                U = results.Up(:, :, snapshot_index);
                V = results.Vp(:, :, snapshot_index);                
                
                %% deform images based on velocity field
                
                % create pixel grid
                [XI,YI] = IMgrid([image_height, image_width], [0 0]);
                XI = XI - 0.5;
                YI = YI - 0.5;
                
                % interpolate displacements onto pixel grid
                UI = VFinterp(X, Y, U, XI, YI, 3);
                VI = VFinterp(X, Y, V, XI, YI, 3);
                
                % translate pixel locations, but
                % since sincBlackmanInterp2 assumes
                % coordinate system is pixel-centered, we need
                % to convert back to index-coordinates for the
                % deform.
                %HOWEVER: we need to use these shifted
                %coordinates later in vector coordinates, so
                %move the -0.5 pixel correction to the call to
                %the sincBlackmanInterp2 function
                XD1 = XI - UI/2;
                YD1 = YI - VI/2;
                XD2 = XI + UI/2;
                YD2 = YI + VI/2;
                                
                im1d = interp2(im1, XD1+0.5, YD1+0.5, 'cubic',0);
                im2d = interp2(im2, XD2+0.5, YD2+0.5, 'cubic',0);

                %% extract interrogation window
                
                % number of vectors along each dimensions
                [grid_rows, grid_cols] = size(X);
                
                % initialize arrays to hold the weights
                w_im_x = nans(grid_rows, grid_cols);
                w_im_y = nans(grid_rows, grid_cols);
                
                w_mc_x = nans(grid_rows, grid_cols);
                w_mc_y = nans(grid_rows, grid_cols);
                
                w_cs_x = nans(grid_rows, grid_cols);
                w_cs_y = nans(grid_rows, grid_cols);
                
                % analysis
                for r = 1:grid_rows
                    for c = 1:grid_cols                        
                        %% identify window extents
                        xmin = X(r, c) - window_size/2 + 1;
                        xmax = X(r, c) + window_size/2;
                        ymin = Y(r, c) - window_size/2 + 1;
                        ymax = Y(r, c) + window_size/2;

                        %% extract image windows
                        
                        %find the image windows
                        zone1 = im1d( max([1 ymin]):min([image_height ymax]),max([1 xmin]):min([image_width xmax]));
                        zone2 = im2d( max([1 ymin]):min([image_height ymax]),max([1 xmin]):min([image_width xmax]));
                        
                        if size(zone1,1)~=window_size || size(zone1,2)~=window_size
                            w1 = zeros(window_size,window_size);
                            w1( 1+max([0 1-ymin]):window_size-max([0 ymax-image_height]), 1+max([0 1-xmin]):window_size-max([0 xmax-image_width]) ) = zone1;
                            zone1 = w1;
                        end
                        if size(zone2,1)~=window_size || size(zone2,2)~=window_size
                            w2 = zeros(window_size,window_size);
                            w2( 1+max([0 1-ymin]):window_size-max([0 ymax-image_height]),1+max([0 1-xmin]):window_size-max([0 xmax-image_width]) ) = zone2;
                            zone2 = w2;
                        end
                                                
                        %% apply windowing function to interrogation window
                        
                        % create window masking filter
                        sfilt1 = windowmask([window_size window_size],[window_resolution window_resolution]);
                        sfilt2 = windowmask([window_size window_size],[window_resolution window_resolution]);
                        
                        % apply the image spatial filter
                        im1_sub = zone1 .* sfilt1;
                        im2_sub = zone2 .* sfilt2;
                        %% identify particles
                        
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
                        [sigma_avg, sigma_std] = calculate_variance_resampling(im1_sub, im2_sub, X(r,c), Y(r,c), U(r,c), V(r,c), x1, y1, x2, y2, d_p*ones(size(x1)), num_particles_remove, num_resampling_trials, uncertainty_flags, intensity_threshold);
                        
                        % temporary fix
                        if sigma_std.mcx > 5 * sigma_std.mcy
                            sigma_std.mcx = sigma_std.mcy;
                        end
                        if sigma_std.csy > 5 * sigma_std.csx
                            sigma_std.csy = sigma_std.csx;
                        end
                        %% calculate weights from bootstrapping results
                        
                                                
                        %% calculate local weights
                        
                        w_im_x(r, c, snapshot_index) = 1/sigma_std.imx^2/(1/sigma_std.imx^2 + 1/sigma_std.mcx^2 + 1/sigma_std.csx^2);
                        w_im_y(r, c, snapshot_index) = 1/sigma_std.imy^2/(1/sigma_std.imy^2 + 1/sigma_std.mcy^2 + 1/sigma_std.csy^2);
                        
                        w_mc_x(r, c, snapshot_index) = 1/sigma_std.mcx^2/(1/sigma_std.imx^2 + 1/sigma_std.mcx^2 + 1/sigma_std.csx^2);
                        w_mc_y(r, c, snapshot_index) = 1/sigma_std.mcy^2/(1/sigma_std.imy^2 + 1/sigma_std.mcy^2 + 1/sigma_std.csy^2);
                        
                        w_cs_x(r, c, snapshot_index) = 1 - (w_im_x(r, c) + w_mc_x(r, c));
                        w_cs_y(r, c, snapshot_index) = 1 - (w_im_y(r, c) + w_mc_y(r, c));
                    end
                end
            end

            %% combine errors and uncertainties
            
            % calculate combined uncertainty
            sigma_U_combined(:, :, snapshot_index) = 1./(w_im_x(:, :, snapshot_index) + w_mc_x(:, :, snapshot_index) + w_cs_x(:, :, snapshot_index)) .* (w_im_x(:, :, snapshot_index) .* results.UIMx(:, :, snapshot_index) + w_mc_x(:, :, snapshot_index) .* results.UMCx(:, :, snapshot_index) + w_cs_x(:, :, snapshot_index) .* sigma_Ud_interp);
            sigma_V_combined(:, :, snapshot_index) = 1./(w_im_y(:, :, snapshot_index) + w_mc_y(:, :, snapshot_index) + w_cs_y(:, :, snapshot_index)) .* (w_im_y(:, :, snapshot_index) .* results.UIMy(:, :, snapshot_index) + w_mc_y(:, :, snapshot_index) .* results.UMCy(:, :, snapshot_index) + w_cs_y(:, :, snapshot_index) .* sigma_Vd_interp);

            % calculate average error
            err_U_combined(:, :, snapshot_index) = 1/2 * (results.err_up(:, :, snapshot_index) + err_Ud_interp);
            err_V_combined(:, :, snapshot_index) = 1/2 * (results.err_vp(:, :, snapshot_index) + err_Vd_interp);

        end

        %% calculate statistics

        % filter out only the valid error estimates
        err_combined_all = [err_U_combined(:); err_V_combined(:)];
        indices = abs(err_combined_all) < max_error_threshold & abs(err_combined_all) > min_error_threshold;
        err_combined_all_valid = nans(size(err_combined_all, 1), size(err_combined_all, 2));
        err_combined_all_valid(indices) = err_combined_all(indices);

        % calculate rms of average error
        err_combined_rms = rms(err_combined_all_valid, 'omitnan');

        % filter out only the valid uncertainty estimates
        sigma_combined_all = [sigma_U_combined(:); sigma_V_combined(:)];
        indices = abs(sigma_combined_all) < max_error_threshold & abs(sigma_combined_all) > min_error_threshold;
        sigma_combined_all_valid = nans(size(sigma_combined_all, 1), size(sigma_combined_all, 2));
        sigma_combined_all_valid(indices) = sigma_combined_all(indices);

        % calculate rms of average uncertainty
        sigma_avg_rms = rms(sigma_combined_all_valid, 'omitnan');

        % calculate pdf of average error
        [N_err_combined, ~] = histcounts(err_combined_all_valid, results_hist.vec, 'normalization', 'pdf');
        N_err_combined = [N_err_combined, 0];

        % calculate pdf of average uncertainty
        [N_sigma_combined, ~] = histcounts(sigma_combined_all_valid, results_hist.vec2, 'normalization', 'pdf');
        N_sigma_combined = [N_sigma_combined, 0];

        %% calculate coverage

        % extract coverage for current case
        switch dataset_name
            case 'PivChal03B'
                results_coverage_current = results_coverage.covg1;        
            case 'PivChal05B'
                results_coverage_current = results_coverage.covg2;
            case 'stagnation_flow'
                results_coverage_current = results_coverage.covg3;
            case 'Vortex_Ring'
                results_coverage_current = results_coverage.covg4;
            case 'Jetdata'
                results_coverage_current = results_coverage.covg5;            
        end

        % calculate combined coverages for MC, IM, and CS
        coverage_mc = 0.5 * (results_coverage_current(1) + results_coverage_current(2));
        coverage_im = 0.5 * (results_coverage_current(3) + results_coverage_current(4));
        coverage_cs = 0.5 * (results_coverage_current(5) + results_coverage_current(6));

        % extract coverages for prana and davis
        coverage_prana = results_coverage_current(7);
        coverage_davis = results_coverage_current(8);

        % ---------------------------------------
        % calculate coverage for combined results
        % ---------------------------------------

        % identify non-nan and valid indices
        non_nan_indices = ~isnan(err_combined_all_valid) & ~isnan(sigma_combined_all_valid);

        % extract non-nan and valid errors and uncertainties
        err_temp = err_combined_all_valid(non_nan_indices);
        sigma_temp = sigma_combined_all_valid(non_nan_indices);

        % calculate coverage (fraction of points with error < uncertainty)
        coverage_combined = sum(abs(err_temp) < abs(sigma_temp))/numel(err_temp) * 100;

        %% calculate coarse histogram for rms error vs uncertainty plot
        
        edges = linspace(0, max_bin, num_bins);
        [~, ~, bin_sigma] = histcounts(sigma_combined_all_valid, edges);

        % calculate binwise rms error and uncertainties
        err_avg_rms_binwise = nans(1, num_bins);
        sigma_avg_rms_binwise = nans(1, num_bins);

        for bin_index = 1:num_bins
            % find uncertainty values lying in the current bin
            sigma_bin_indices = find(bin_sigma == bin_index);
            % calculate rms of the uncertainties of the values in the
            % current bin
            sigma_avg_rms_binwise(bin_index) = rms(sigma_combined_all_valid(sigma_bin_indices), 'omitnan');
            % calculate rms of the errors of the measurements in the
            % current bin
            err_avg_rms_binwise(bin_index) = rms(err_combined_all_valid(sigma_bin_indices), 'omitnan');
        end

        
        %% display results

        % =================================
        % Prana, IM and MC contours
        % =================================

        figure
        subplot(3,1,1)
        contourf(results.Xp(:, :, 1), results.Yp(:, :, 1), sqrt(results.Up(:, :, 1).^2 + results.Vp(:, :, 1).^2), displacement_contour_levels, 'edgecolor', 'none');
        colormap(flipud(gray))
        annotate_image(gcf, gca);
        h2 = colorbar;
        caxis([displacement_color_min displacement_color_max])
        title(h2, '\Delta x (pix.)')
        title('Displacement')

        subplot(3,1,2)
        contourf(results.Xp(:, :, 1), results.Yp(:, :, 1), sqrt(results.UIMx(:, :, 1).^2 + results.UIMy(:, :, 1).^2), uncertainty_contour_levels, 'edgecolor', 'none');
        colormap(flipud(gray))
        annotate_image(gcf, gca);
        h2 = colorbar;
        caxis([uncertainty_color_min uncertainty_color_max])
        title(h2, '\sigma_{\Delta x} (pix.)')
        title({'Image'; 'Matching'})

        subplot(3,1,3)
        contourf(results.Xp(:, :, 1), results.Yp(:, :, 1), sqrt(results.UMCx(:, :, 1).^2 + results.UMCy(:, :, 1).^2), uncertainty_contour_levels, 'edgecolor', 'none');
        colormap(flipud(gray))
        annotate_image(gcf, gca);
        h2 = colorbar;
        caxis([uncertainty_color_min uncertainty_color_max])
        title(h2, '\sigma_{\Delta x} (pix.)')
        title({'Moment of'; 'Correlation'})
        set(gcf, 'Position', [132    34   646   851]);

        if save_figures
            save_figure_to_png_eps_fig(figure_write_directory, 'prana-im-mc-contours', [true, false, false]);    
        end

        % =================================
        % Davis and CS contours
        % =================================

        figure
        subplot(2, 1, 1)
        contourf(results.Xd(:, :, 1), results.Yd(:, :, 1), sqrt(results.Ud(:, :, 1).^2 + results.Vd(:, :, 1).^2), displacement_contour_levels, 'edgecolor', 'none');
        colormap(flipud(gray))
        annotate_image(gcf, gca);
        h2 = colorbar;
        caxis([displacement_color_min displacement_color_max])
        title(h2, '\Delta x (pix.)')
        title('Displacement')

        subplot(2, 1, 2)
        contourf(results.Xd(:, :, 1), results.Yd(:, :, 1), sqrt(results.UCSx(:, :, 1).^2 + results.UCSy(:, :, 1).^2), uncertainty_contour_levels, 'edgecolor', 'none');
        colormap(flipud(gray))
        annotate_image(gcf, gca);
        h2 = colorbar;
        caxis([uncertainty_color_min uncertainty_color_max])
        title(h2, '\sigma_{\Delta x} (pix.)')
        title({'Correlation'; 'Statistics'})

        set(gcf, 'Position', [132         310        646         575]);

        if save_figures
            save_figure_to_png_eps_fig(figure_write_directory, 'davis-cs-contours', [true, false, false]);    
        end

        % =================================
        % Davis and CS contours Interpolated
        % =================================

        figure
        subplot(2, 2, 1)
        contourf(results.Xd(:, :, 1), results.Yd(:, :, 1), sqrt(results.Ud(:, :, end).^2 + results.Vd(:, :, end).^2), displacement_contour_levels, 'edgecolor', 'none');
        colormap(flipud(gray))
        annotate_image(gcf, gca);
        h2 = colorbar;
        caxis([displacement_color_min displacement_color_max])
        title(h2, '\Delta x (pix.)')
        title('Displacement')

        subplot(2, 2, 2)
        contourf(results.Xd(:, :, 1), results.Yd(:, :, 1), sqrt(results.UCSx(:, :, end).^2 + results.UCSy(:, :, end).^2), uncertainty_contour_levels, 'edgecolor', 'none');
        colormap(flipud(gray))
        annotate_image(gcf, gca);
        h2 = colorbar;
        caxis([uncertainty_color_min uncertainty_color_max])
        title(h2, '\sigma_{\Delta x} (pix.)')
        title({'Correlation'; 'Statistics'})


        subplot(2, 2, 3)
        contourf(results.Xp(:, :, 1), results.Yp(:, :, 1), sqrt(Ud_interp.^2 + Vd_interp.^2), displacement_contour_levels, 'edgecolor', 'none');
        colormap(flipud(gray))
        annotate_image(gcf, gca);
        h2 = colorbar;
        caxis([displacement_color_min displacement_color_max])
        title(h2, '\Delta x (pix.)')
        title('Displacement, Interpolated')

        subplot(2, 2, 4)
        contourf(results.Xp(:, :, 1), results.Yp(:, :, 1), sqrt(sigma_Ud_interp.^2 + sigma_Vd_interp.^2), uncertainty_contour_levels, 'edgecolor', 'none');
        colormap(flipud(gray))
        annotate_image(gcf, gca);
        h2 = colorbar;
        caxis([uncertainty_color_min uncertainty_color_max])
        title(h2, '\sigma_{\Delta x} (pix.)')
        title({'Correlation'; 'Statistics, Interp.'})

        % set(gcf, 'Position', [132         310        646         575]);
        set(gcf, 'Position', [51         367        1217         549]);

        if save_figures
            save_figure_to_png_eps_fig(figure_write_directory, 'davis-cs-interp-contours', [true, false, false]);    
        end

        % =================================
        % Combined uncertainty contours
        % =================================

        figure
        subplot(2, 1, 1)
        contourf(results.Xp(:, :, 1), results.Yp(:, :, 1), sqrt(U_avg(:, :, 1).^2 + V_avg(:, :, 1).^2), displacement_contour_levels, 'edgecolor', 'none');
        colormap(flipud(gray))
        annotate_image(gcf, gca);
        h2 = colorbar;
        caxis([displacement_color_min displacement_color_max])
        title(h2, '\Delta x (pix.)')
        title('Displacement, Averaged')

        subplot(2,1,2)
        contourf(results.Xp(:, :, 1), results.Yp(:, :, 1), sqrt(sigma_U_combined(:, :, 1).^2 + sigma_V_combined(:, :, 1).^2), uncertainty_contour_levels, 'edgecolor', 'none');
        colormap(flipud(gray))
        annotate_image(gcf, gca);
        h2 = colorbar;
        caxis([uncertainty_color_min uncertainty_color_max])
        title(h2, '\sigma_{\Delta x} (pix.)')
        title({'Uncertainty'; 'Combined'})

        set(gcf, 'Position', [132         310        646         575]);

        if save_figures
            save_figure_to_png_eps_fig(figure_write_directory, 'combined-uncertainty-contours', [true, false, false]);    
        end

        close all;
        %% aggregate results for this uncertainty combination scheme
        
        % allocate memory
        if uncertainty_combination_method_index == 1
            N_sigma_combined_all_schemes = nans(numel(N_sigma_combined), num_uncertainty_combination_methods);
            sigma_avg_rms_all_schemes = nans(numel(sigma_avg_rms), num_uncertainty_combination_methods);
            err_avg_rms_binwise_all_schemes = nans(numel(err_avg_rms_binwise), num_uncertainty_combination_methods);
            sigma_avg_rms_binwise_all_schemes = nans(numel(sigma_avg_rms_binwise), num_uncertainty_combination_methods);
            coverage_combined_all_schemes = nans(numel(coverage_combined), num_uncertainty_combination_methods);
        end
        
        N_sigma_combined_all_schemes(:, uncertainty_combination_method_index) = N_sigma_combined';
        sigma_avg_rms_all_schemes(uncertainty_combination_method_index) = sigma_avg_rms;
        err_avg_rms_binwise_all_schemes(:, uncertainty_combination_method_index) = err_avg_rms_binwise';
        sigma_avg_rms_binwise_all_schemes(:, uncertainty_combination_method_index) = sigma_avg_rms_binwise';
        coverage_combined_all_schemes(uncertainty_combination_method_index) = coverage_combined;
    end
    
    %% plot results
    
    % =================================
    % Error and Uncertainty histograms
    % =================================

    figure
    % prana error
    l1 = plot(results_hist.sigprana, results_hist.ll, 'k--');
    hold on
    % davis error
    l2 = plot(results_hist.sigdavis, results_hist.ll, 'm--');
    % average error
    l3 = plot(err_combined_rms*[1 1], [0 30], 'g--');
    % MC
    l4 = plot(results_hist.vec2, results_hist.Nmc, 'r');
    plot(results_hist.sigMC, results_hist.ll, 'r')
    % IM
    l5 = plot(results_hist.vec2, results_hist.Nim, 'c');
    plot(results_hist.sigIM, results_hist.ll, 'c')
    % CS
    l6 = plot(results_hist.vec2, results_hist.Ncs, 'm');
    plot(results_hist.sigCS, results_hist.ll, 'm')
    
    l_wt_none = plot(results_hist.vec2, N_sigma_combined_all_schemes(:, 1), 'g');
    plot(sigma_avg_rms_all_schemes(:, 1) * ones(1,100), linspace(0, y_max, 100), 'g')
    
    l_wt_global = plot(results_hist.vec2, N_sigma_combined_all_schemes(:, 2), 'g--');
    plot(sigma_avg_rms_all_schemes(:, 2) * ones(1,100), linspace(0, y_max, 100), 'g--')
    
    l_wt_local = plot(results_hist.vec2, N_sigma_combined_all_schemes(:, 3), 'g:');
    plot(sigma_avg_rms_all_schemes(:, 3) * ones(1,100), linspace(0, y_max, 100), 'g:')
        
    legend([l_im, l_mc, l_cs, l_wt_none, l_wt_global, l_wt_local, l_true], ...
        {'\sigma_{IM}', '\sigma_{MC}', '\sigma_{CS}', 'Unweighted', 'Global', 'Local', 'RMS, True'}, 'fontsize', 12)

%     % average uncertainty
%     l7 = plot(results_hist.vec2, N_sigma_combined, 'g');
%     plot(sigma_avg_rms*[1 1], [0 30], 'g')

    ylim([0 100])
    xlabel('(pix.)')
    ylabel('Count')
%     legend([l1, l2, l3, l4, l5, l6, l7], {'e_{prana}', 'e_{davis}', 'e_{comb}', '\sigma_{MC}', '\sigma_{IM}', '\sigma_{CS}', '\sigma_{comb}'}, 'location', 'northeast')
    legend([l1, l2, l3, l4, l5, l6, l_wt_none, l_wt_global, l_wt_local], {'e_{prana}', 'e_{davis}', 'e_{comb}', '\sigma_{MC}', '\sigma_{IM}', '\sigma_{CS}', 'Unweighted', 'Global', 'Local'}, 'location', 'northeast')
    return;
    set(gcf, 'Position', [336   573   429   321])
    
    if save_figures
        save_figure_to_png_eps_fig(figure_write_directory, 'histograms-all-schemes', [true, false, false]);    
    end

    return;
    %%

    % ==============================
    % RMS Error vs RMS Uncertainty
    % ==============================

    figure
    % Error
    plot([0 max_bin], [0 max_bin], 'k--');
    hold on
    % MC
    l1 = plot(results_rms.rmserru1, results_rms.rmsMC, 'ro-', 'markerfacecolor', 'r');
    % IM
    l2 = plot(results_rms.rmserru2, results_rms.rmsIM, 'co-', 'markerfacecolor', 'c');
    % CS
    l3 = plot(results_rms.rmserru3, results_rms.rmsCS, 'mo-', 'markerfacecolor', 'm');
    % avg
    l4 = plot(err_avg_rms_binwise, sigma_avg_rms_binwise, 'go-', 'markerfacecolor', 'g');

    grid on
    axis equal
    axis([0 max_bin 0 max_bin])
    xlabel('RMS Error (pix.)')
    ylabel('RMS Uncertainty (pix.)')
    legend([l1, l2, l3, l4], {'MC', 'IM', 'CS', 'Combined'}, 'location', 'northwest')
    set(gcf, 'Position', [336   511   488   383])

    if save_figures
        save_figure_to_png_eps_fig(figure_write_directory, 'rms-error-uncertainty', [true, false, false]);    
    end

    %% 

    % ====================================
    % normalized weights across snapshots
    % ====================================

    figure
    % x
    subplot(1, 2, 1)
    plot(w_mc_x./(w_im_x + w_mc_x + w_cs_x), 'ro', 'markerfacecolor', 'r')
    hold on
    plot(w_im_x./(w_im_x + w_mc_x + w_cs_x), 'co', 'markerfacecolor', 'c')
    plot(w_cs_x./(w_im_x + w_mc_x + w_cs_x), 'mo', 'markerfacecolor', 'm')
    ylim([0 1])
    xlabel('Snapshot #')
    ylabel('Weight (x)')

    legend('MC', 'IM', 'CS', 'location', 'north', 'orientation', 'horizontal')
    legend boxon

    % y
    subplot(1, 2, 2)
    plot(w_mc_y./(w_im_y + w_mc_y + w_cs_y), 'ro', 'markerfacecolor', 'r')
    hold on
    plot(w_im_y./(w_im_y + w_mc_y + w_cs_y), 'co', 'markerfacecolor', 'c')
    plot(w_cs_y./(w_im_y + w_mc_y + w_cs_y), 'mo', 'markerfacecolor', 'm')
    ylim([0 1])
    xlabel('Snapshot #')
    ylabel('Weight (y)')

    legend('MC', 'IM', 'CS', 'location', 'north', 'orientation', 'horizontal')
    legend boxon
    set(gcf, 'Position', [147   448   779   386])

    if save_figures
        save_figure_to_png_eps_fig(figure_write_directory, 'normalized weights', [true, false, false]);    
    end

    %%

    % ====================================
    % coverage
    % ====================================

    figure

    X = categorical({'MC','IM','CS','Comb'});
    X = reordercats(X,{'MC','IM','CS','Comb'});

    Y = [coverage_mc, coverage_im, coverage_cs, coverage_combined];

    bar(X,Y); %,'FaceColor','flat');
    hold on
    l1 = plot(X, coverage_prana * [1, 1, 1, 1], 'r');
    l2 = plot(X, coverage_davis * [1, 1, 1, 1], 'g');
    ylim([0 100])
    ylabel('Coverage (%)')
    legend([l1, l2], {'Target, Prana'; 'Target, Davis'}, 'location', 'northoutside', 'Orientation', 'horizontal')
    set(gcf, 'Position', [336   573   429   321])

    if save_figures
        save_figure_to_png_eps_fig(figure_write_directory, 'coverage', [true, false, false]);    
    end

    %% close all figures
    close all;
    
end
clear
close all
clc


restoredefaultpath;
addpath(genpath('/scratch/shannon/c/aether/Projects/BOS/general-codes/matlab-codes/'));
setup_default_settings;
% addpath ../prana-uncertainty-average-dc-new-im-2-save-deform-cs/
addpath ../prana/
addpath ../general-codes/

% ========================
%% processing settings
% ========================
% directory containing images
top_image_directory = fullfile('/scratch/shannon/a/aether/Projects/PIV_PTV_Uncertainty_Quantification/Planar_uncertainty_work/Images/experiment/');

% array of data set names
dataset_name_array = {'PivChal03B'; 'PivChal05B'; 'stagnation_flow'; 'Vortex_Ring'; 'Jetdata'};
num_datasets = numel(dataset_name_array);

% directory containing sample job files
sample_job_file_directory = '/scratch/shannon/a/aether/Projects/PIV_PTV_Uncertainty_Quantification/uncertainty-combination/analysis/data/sample-job-files/';

% directory where results of this analysis are to be saved
% top_write_directory = fullfile('/scratch/shannon/a/aether/Projects/PIV_PTV_Uncertainty_Quantification/uncertainty-combination/analysis/results/planar-dataset/', 'experiment-new-grid-buffer');
top_write_directory = fullfile('/scratch/shannon/a/aether/Projects/PIV_PTV_Uncertainty_Quantification/uncertainty-combination/analysis/results/planar-dataset/', 'experiment-new');
if ~exist(top_write_directory, 'dir')
    mkdir(top_write_directory);
end

% array of uncertainty methods
individual_method_array = {'IM'; 'MC'; 'CS'};

% pass number
pass_number = 4;

% number of trials for comparing
num_trials = 1e2;

% prana processing?
process_prana = 0;

% save figures?
save_figures = 1;

% ========================
%% run processing
% ========================
for window_size_index = 1:2
    fprintf('WS: %d\n', window_size_index);        
    for dataset_index = 1:num_datasets
        % extract data set name
        dataset_name = dataset_name_array{dataset_index};
        fprintf('dataset name: %s\n', dataset_name);
        
        % image directory for current data set
        current_image_directory = fullfile(top_image_directory, dataset_name);
        
        % results directory for current data set
        current_results_directory = fullfile(top_write_directory, dataset_name, ['WS' num2str(window_size_index)]);
        mkdir_c(current_results_directory);
        
        % ====================================
        % original processing
        % ====================================
        fprintf('Original Processing\n');
        % load job file
        sample_job_file = load(fullfile(current_results_directory, 'jobfile.mat'));
        
        % extract data structure containing job parameters
        Data = sample_job_file.Data;

        % directory to save vectors
        Data.outdirec = fullfile(current_results_directory, 'vectors-temp');
        mkdir_c(Data.outdirec);

        % number of files to be processed
        Data.imfend = '1';
                
        % turn off parallel processing if just one image pair needs to be
        % processed
        if str2double(Data.imfend) == 1
            Data.par = '0';
        end
        
        % call prana
        if process_prana
            pranaPIVcode(Data);
        end

        % ====================================
        % new processing
        % ====================================
        fprintf('New Processing\n');

        snapshot_index = 1;

        % get results listing
        [results_all, num_results] = load_directory_data(Data.outdirec, ['PIV_pass' num2str(pass_number) '*.mat']);
        % load results
        % results = load(fullfile(Data.outdirec, 'PIV_pass4_0001.mat'));
        results = results_all(snapshot_index);

        % extract current pass settings
        current_pass_settings = Data.(['PIV' num2str(pass_number)]);

        % extract size of the co-ordinates
        [num_rows, num_cols] = size(results.X);
        
        % ------------------------
        %% load listing of deformed images
        % ------------------------
        % directory containing deformed images for current data set
        deformed_images_directory = fullfile(Data.outdirec, 'imDeform');
        
        % get list of im1 files in the directory
        files_im1 = get_directory_listing(deformed_images_directory, 'PIV*im1d*.mat');
        % get list of im2 files in the directory
        files_im2 = get_directory_listing(deformed_images_directory, 'PIV*im2d*.mat');

        % load deformed images for the first frame
        % im1 = imread(fullfile(files_im1(1).folder, files_im1(1).name));
        % im1d = double(flipud(im1));
        saved_img = load(fullfile(files_im1(snapshot_index).folder, files_im1(snapshot_index).name));
        im1d = saved_img.im1d;

        % load deformed images for the second frame
        % im2 = imread(fullfile(files_im2(1).folder, files_im2(1).name));
        % im2d = double(flipud(im2));
        saved_img = load(fullfile(files_im2(snapshot_index).folder, files_im2(snapshot_index).name));        
        im2d = saved_img.im2d;

        % extract image height and width
        [image_height, image_width] = size(im1d);

        % extract window size
        [window_size_x, window_size_y] = extract_window_size(current_pass_settings);

        % ------------------------
        % loop through grid points
        % ------------------------
        r_trials = nans(1, num_trials);
        c_trials = nans(1, num_trials);
        imx = nans(2, num_trials);
        imy = nans(2, num_trials);
        mcx = nans(2, num_trials);
        mcy = nans(2, num_trials);
        csx = nans(2, num_trials);
        csy = nans(2, num_trials);

        rng(0);
        for trial_index = 1:num_trials
            fprintf('trial_index: %d\n', trial_index);
            % % extract individual uncertainties for this grid point
            % unc_original = extract_planar_uncertainties(results.uncertainty2D, r, c);
            
            % calculate row and column indicies
            r = randi(num_rows);
            c = randi(num_cols);

            % store values
            r_trials(trial_index) = r;
            c_trials(trial_index) = c;

            % extract current grid point co-ordinates
            X = results.X(r, c);
            Y = results.Y(r, c);

            % extract original uncertainties at these points
            imx(1, trial_index) = results.uncertainty2D.Uimx(r, c);
            imy(1, trial_index) = results.uncertainty2D.Uimy(r, c);
            mcx(1, trial_index) = results.uncertainty2D.MCx(r, c);
            mcy(1, trial_index) = results.uncertainty2D.MCy(r, c);
            csx(1, trial_index) = results.uncertainty2D.Ucsx(r, c);
            csy(1, trial_index) = results.uncertainty2D.Ucsy(r, c);
            
            % extract interrogation window
            im1_sub = extract_interrogation_window(im1d, X, Y, [window_size_x, window_size_y]);
            im2_sub = extract_interrogation_window(im2d, X, Y, [window_size_x, window_size_y]);
            
            % calculate planar uncertainty on the original interrogation windows
            % [unc_sub_current, snr_metric_sub_current] = calculate_planar_uncertainties(im1_sub, im2_sub, [window_size_x, window_size_y], [window_resolution_x, window_resolution_y], individual_method_array, uncertainty_flags);
            [unc_sub_current, snr_metric_sub_current, U_sub_current, V_sub_current] = calculate_planar_uncertainties(im1_sub, im2_sub, current_pass_settings, individual_method_array);

            % extract uncertainties
            imx(2, trial_index) = unc_sub_current{1}.x;
            imy(2, trial_index) = unc_sub_current{1}.y;

            mcx(2, trial_index) = unc_sub_current{2}.x;
            mcy(2, trial_index) = unc_sub_current{2}.y;
            
            csx(2, trial_index) = unc_sub_current{3}.x;
            csy(2, trial_index) = unc_sub_current{3}.y;
        end

        % save results
        % save(fullfile(Data.outdirec, 'new-uncertainties.mat'), 'imx', 'imy', 'mcx', 'mcy', 'csx', 'csy', 'r_trials', 'c_trials');

        % ====================================
        % plot comparison
        % ====================================
        plot_original_new_uncertainties_02(imx, imy, mcx, mcy, csx, csy);
        return;
        if save_figures
            save_figure_to_png_svg_fig(Data.outdirec, 'uncertainty-comparison-original-new', [1, 0, 0]);
        end
    end
end


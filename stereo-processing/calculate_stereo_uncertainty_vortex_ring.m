%Uncertainty Propagation for Stereo Geometric Reconstruction(Willert
%Method)for vortex ring, where planar uncertainties are evaluated using CS
%method in Davis
clc;
clear;
close all;

restoredefaultpath;
addpath(genpath('/scratch/shannon/c/aether/Projects/BOS/general-codes/matlab-codes/'));
addpath('../prana-uncertainty-average-dc-new-im-2-save-deform-cs/');
addpath('../stereo_uncertainty_codes_packaged/');
setup_default_settings;

% ============================
% experiment settings
% ============================
% seconds per frame
spf = 0.001;
% No. of frames
num_frames = 50;   
% Starting frame
fstart = 24;
% camera numbers
camera_numbers = [2, 4];
num_cameras = numel(camera_numbers);

% ============================
% directory settings
% ============================
% top level read directory
top_read_directory = '/scratch/shannon/a/aether/Projects/PIV_PTV_Uncertainty_Quantification/Stereo_uncertainty_work/';
% top level write directory
% top_write_directory = '/scratch/shannon/a/aether/Projects/PIV_PTV_Uncertainty_Quantification/uncertainty-combination/analysis/results/stereo-dataset/cam13/With_selfcal/';
top_write_directory = fullfile('/scratch/shannon/a/aether/Projects/PIV_PTV_Uncertainty_Quantification/uncertainty-combination/analysis/results/stereo-dataset/', ...
                                ['cam' num2str(camera_numbers(1)) num2str(camera_numbers(2))], '/With_selfcal/');

% directory containing results
code_package_directory = '../stereo_uncertainty_codes_packaged/';
% 2d uncertainty methods
individual_methods = {'IM'; 'MC'; 'CS'};
num_individual_methods = numel(individual_methods);

% ============================
% processing settings
% ============================
% stereo reconstruction type
rectype = 'Willert';
% pass number of results
pass_index = 4;

% ============================
% plot settings
% ============================
% save figure? (true/false)
save_figures = 1;
% directory to save figures
figure_save_directory = fullfile(top_write_directory, 'figures', ['pass' num2str(pass_index)]);
mkdir_c(figure_save_directory);
% uncertainty limits
uncertainty_limits = [0, 1];
% contour levels for uncertainty
uncertainty_contour_levels = linspace(uncertainty_limits(1), uncertainty_limits(2), 100);

% ============================
% Load Stereo calibration job
% ============================
fprintf('loading stereo calibration job...\n');
stereojob = load(fullfile(code_package_directory,'Prana_stereo_job_after_selfcal.mat'));
stdjob = stereojob.stdjob;

% extract Calibration job
caljob = stdjob.caljobfile.datasave.caljob;
calmat = [caljob.aXcam1 caljob.aYcam1 caljob.aXcam2 caljob.aYcam2];

% extract 2D prana job
planarjob = stdjob.pranajob.Data;

% Load disparity field between two camera using ensemble correlation
dispfield = load(fullfile(code_package_directory, 'VR13_scc_pass2_00025.mat'));

% ============================
% load planar fields
% ============================
fprintf('loading planar results...\n');
fprintf('cam1\n');
% results directory for cam 1
cam1_results_directory = fullfile(top_write_directory, rectype, ['Camera', num2str(caljob.camnumber(1)), filesep]);
% load mat files
[cam1_results, num_cam1_results] = load_directory_data(fullfile(cam1_results_directory, 'vectors'), ['VR*pass' num2str(pass_index, '%d') '*.mat']);
fprintf('cam2\n');
% results directory for cam 2
cam2_results_directory = fullfile(top_write_directory, rectype, ['Camera', num2str(caljob.camnumber(2)), filesep]);
% load mat files
[cam2_results, num_cam2_results] = load_directory_data(fullfile(cam2_results_directory, 'vectors'), ['VR*pass' num2str(pass_index, '%d') '*.mat']);

% ============================
% load reconstructed fields
% ============================
fprintf('loading stereo results...\n');
% directory containing processed result
stereo_results_directory = fullfile(top_write_directory, rectype, 'Camera1Camera2_3Cfields',filesep);
% load mat files
[stereo_vectors_results, num_stereo_vectors_results] = load_directory_data(fullfile(stereo_results_directory, 'vectors'), ['piv*pass_' num2str(pass_index, '%d') '*.mat']);

% Load Davis individual camera fields and stereo reconstructed fields 
Davissol = load(fullfile(code_package_directory, 'Davis_processed_result.mat'));
Davis13 = Davissol.Davis13;

% ============================
% Load true solution
% ============================
% % true_solution = load(fullfile())
% true_solution = load(fullfile(code_package_directory, 'Vortex_Ring_True_Solution_interp_on_solution_grid.mat'));
% % extract true solution velocities
% Ut = true_solution.Ut;
% Vt = true_solution.Vt;
% Wt = true_solution.Wt;
true_solution  = load_directory_data('/scratch/shannon/a/aether/Projects/PIV_PTV_Uncertainty_Quantification/Stereo_uncertainty_work/Vortex_ring_true_solution/Eval_1/', 'E_eval1_TrueSolution*'); 
[xt, yt] = meshgrid(-34:0.45:20, -27:0.45:27);

% ============================
% Extract calibration details
% ============================
% Magnification in mm/pix
scaling = stdjob.scaling.wil;
mx = scaling.xscale;
my = scaling.yscale;
% Here the results are for camera1camera3 pair
cam1 = '1';
cam2 = '3';
if strcmp(cam1,'1') && strcmp(cam2,'3')
    mz = scaling.yscale;
elseif strcmp(cam1,'2') && strcmp(cam2,'4')
    mz = scaling.xscale;
end

% ============================
%% Find uncertainty in angles 
% ============================
fprintf('calculating uncertainty in angles\n');
% directory which has uncertainty in disparity field
disp_uncertainty_filedir = code_package_directory;
% planarjob.outdirec = '/scratch/shannon/a/aether/Projects/PIV_PTV_Uncertainty_Quantification/Stereo_uncertainty_work/Results/Prana_process/results/cam13/With_selfcal/';
planarjob.outdirec = top_write_directory;

[Un_alpha1,Un_alpha2,Un_beta1,Un_beta2,tanalpha1,tanalpha2,tanbeta1,tanbeta2] = stereo_angle_uncertainty(caljob,planarjob,dispfield,disp_uncertainty_filedir);

% ============================
%% calculate stereo uncertainty
% ============================
fprintf('Begin uncertainty propagation\n');
% Initialization
[Sx,Sy] = size(tanalpha1);
Us = zeros(Sx,Sy,num_frames); Vs = zeros(Sx,Sy,num_frames); Ws = zeros(Sx,Sy,num_frames);
Ut = zeros(Sx,Sy,num_frames); Vt = zeros(Sx,Sy,num_frames); Wt = zeros(Sx,Sy,num_frames);
err_u = zeros(Sx,Sy,num_frames); err_v = zeros(Sx,Sy,num_frames); err_w = zeros(Sx,Sy,num_frames);
unc1_u = {zeros(Sx,Sy,num_frames); zeros(Sx,Sy,num_frames); zeros(Sx,Sy,num_frames)};
unc1_v = {zeros(Sx,Sy,num_frames); zeros(Sx,Sy,num_frames); zeros(Sx,Sy,num_frames)};
unc2_u = {zeros(Sx,Sy,num_frames); zeros(Sx,Sy,num_frames); zeros(Sx,Sy,num_frames)};
unc2_v = {zeros(Sx,Sy,num_frames); zeros(Sx,Sy,num_frames); zeros(Sx,Sy,num_frames)};
unc_st_u = {zeros(Sx,Sy,num_frames); zeros(Sx,Sy,num_frames); zeros(Sx,Sy,num_frames)};
unc_st_v = {zeros(Sx,Sy,num_frames); zeros(Sx,Sy,num_frames); zeros(Sx,Sy,num_frames)};
unc_st_w = {zeros(Sx,Sy,num_frames); zeros(Sx,Sy,num_frames); zeros(Sx,Sy,num_frames)};

% loop through time steps
for frame_index=1:num_frames
    fprintf('frame_index = %d of %d\n', frame_index, num_frames);
    % ============================
    %% Load reconstructed solution
    % ============================
    % convert velocities to pixel units
    % Us(:,:,frame_index)=(spf/mx)*1e3*Davis13.Uw(1:end-1,1:end-1,frame_index);
    % Vs(:,:,frame_index)=(spf/my)*1e3*Davis13.Vw(1:end-1,1:end-1,frame_index);
    % Ws(:,:,frame_index)=(spf/mz)*1e3*Davis13.Ww(1:end-1,1:end-1,frame_index);
    Us(:,:,frame_index)=(spf/mx)*1e3*stereo_vectors_results{frame_index}.U;
    Vs(:,:,frame_index)=(spf/mx)*1e3*stereo_vectors_results{frame_index}.V;
    Ws(:,:,frame_index)=(spf/mx)*1e3*stereo_vectors_results{frame_index}.W;
    
    % convert positions to physical units
    % xgrid=1000*Davis13.Xw(1:end-1,1:end-1);
    % ygrid=1000*Davis13.Yw(1:end-1,1:end-1);
    % zgrid=1000*Davis13.Zw(1:end-1,1:end-1);
    xgrid=1000*stereo_vectors_results{frame_index}.X;
    ygrid=1000*stereo_vectors_results{frame_index}.Y;
    zgrid=1000*stereo_vectors_results{frame_index}.Z;

    % ============================
    %% Get error fields
    % ============================
    % interpolate true solution on to measurement grid
    ut = spf/mx * 1000 * true_solution{frame_index}.velocity_data.x_velocity;
    vt = spf/mx * 1000 * true_solution{frame_index}.velocity_data.y_velocity;
    wt = spf/mx * 1000 * true_solution{frame_index}.velocity_data.z_velocity;
    Ut(:, :, frame_index) = interp2(xt, yt, ut, xgrid, ygrid, 'linear', 0);
    Vt(:, :, frame_index) = interp2(xt, yt, vt, xgrid, ygrid, 'linear', 0);
    Wt(:, :, frame_index) = interp2(xt, yt, wt, xgrid, ygrid, 'linear', 0);
    
    % Ut_interp = interp2(true_solution.xgrid, true_solution.ygrid, Ut(:,:,frame_index), xgrid, ygrid);
    % Vt_interp = interp2(true_solution.xgrid, true_solution.ygrid, Vt(:,:,frame_index), xgrid, ygrid);
    % Wt_interp = interp2(true_solution.xgrid, true_solution.ygrid, Wt(:,:,frame_index), xgrid, ygrid);
    err_u(:,:,frame_index)=(Ut(:,:,frame_index) - Us(:,:,frame_index));
    err_v(:,:,frame_index)=(Vt(:,:,frame_index) - Vs(:,:,frame_index));
    err_w(:,:,frame_index)=(Wt(:,:,frame_index) - Ws(:,:,frame_index));

    % ============================
    %% Load cam1cam2 2d fields
    % ============================
    cam1_result_current = cam1_results{frame_index};
    cam2_result_current = cam2_results{frame_index};

    % get 2d velocities
    % U1 = Davis13.Us1(1:end-1,1:end-1,frame_index);
    % V1 = Davis13.Vs1(1:end-1,1:end-1,frame_index);
    % X1 = Davis13.X(1:end-1,1:end-1);
    % Y1 = Davis13.Y(1:end-1,1:end-1);
    U1 = cam1_result_current.U;
    V1 = cam1_result_current.V;
    X1 = cam1_result_current.X;
    Y1 = cam1_result_current.Y;

    % U2 = Davis13.Us2(1:end-1,1:end-1,frame_index);
    % V2 = Davis13.Vs2(1:end-1,1:end-1,frame_index);
    % X2 = Davis13.X(1:end-1,1:end-1);
    % Y2 = Davis13.Y(1:end-1,1:end-1);
    U2 = cam2_result_current.U;
    V2 = cam2_result_current.V;
    X2 = cam2_result_current.X;
    Y2 = cam2_result_current.Y;
    
    for method_index = 1:num_individual_methods
        % ============================
        % Extract 2d uncertainties
        % ============================
        % Camera1
        % unc1_u(:,:,frame_index) = Davis13.UCSx1(1:end-1,1:end-1,frame_index);
        % unc1_v(:,:,frame_index) = Davis13.UCSy1(1:end-1,1:end-1,frame_index);
        % Camera2
        % unc2_u(:,:,frame_index) = Davis13.UCSx2(1:end-1,1:end-1,frame_index);
        % unc2_v(:,:,frame_index) = Davis13.UCSy2(1:end-1,1:end-1,frame_index);
        if method_index == 1
            % IM
            unc1_u{method_index}(:, :, frame_index) = cam1_result_current.uncertainty2D.Uimx;
            unc1_v{method_index}(:, :, frame_index) = cam1_result_current.uncertainty2D.Uimy;
            unc2_u{method_index}(:, :, frame_index) = cam2_result_current.uncertainty2D.Uimx;
            unc2_v{method_index}(:, :, frame_index) = cam2_result_current.uncertainty2D.Uimy;
        elseif method_index == 2
            % MC
            unc1_u{method_index}(:, :, frame_index) = sqrt(cam1_result_current.uncertainty2D.MCx.^2); % - cam1_results(frame_index).uncertainty2D.biasx.^2);
            unc1_v{method_index}(:, :, frame_index) = sqrt(cam1_result_current.uncertainty2D.MCy.^2); % - cam1_results(frame_index).uncertainty2D.biasy.^2);
            unc2_u{method_index}(:, :, frame_index) = sqrt(cam2_result_current.uncertainty2D.MCx.^2); % - cam2_results(frame_index).uncertainty2D.biasx.^2);
            unc2_v{method_index}(:, :, frame_index) = sqrt(cam2_result_current.uncertainty2D.MCy.^2); % - cam2_results(frame_index).uncertainty2D.biasy.^2);            
            % % bias correction
            % unc1_u{method_index}(:, :, frame_index) = cam1_results(frame_index).uncertainty2D
        else
            % CS
            unc1_u{method_index}(:, :, frame_index) = cam1_result_current.uncertainty2D.Ucsx;
            unc1_v{method_index}(:, :, frame_index) = cam1_result_current.uncertainty2D.Ucsy;
            unc2_u{method_index}(:, :, frame_index) = cam2_result_current.uncertainty2D.Ucsx;
            unc2_v{method_index}(:, :, frame_index) = cam2_result_current.uncertainty2D.Ucsy;
        end
        
        % ============================
        %% Stereo Uncertainty Propagation
        % ============================
        % extract uncertainties
        Unu1 = unc1_u{method_index}(:,:,frame_index);
        Unu2 = unc2_u{method_index}(:,:,frame_index);
        Unv1 = unc1_v{method_index}(:,:,frame_index);
        Unv2 = unc2_v{method_index}(:,:,frame_index);

        % propagate uncertainties
        [Un_u1, Un_v1, Un_w1, JU, JV, JW] = stereo_uncertainty_propagation(Unu1, Unv1, Unu2, Unv2, U1, V1, U2, V2, Ws(:,:,frame_index), ...
                                            Un_alpha1, Un_alpha2, Un_beta1, Un_beta2, tanalpha1, tanalpha2, tanbeta1, tanbeta2, mx, my);
        
        % save propagated uncertainties
        unc_st_u{method_index}(:,:,frame_index) = Un_u1;
        unc_st_v{method_index}(:,:,frame_index) = Un_v1;
        unc_st_w{method_index}(:,:,frame_index) = Un_w1;
    end
end

% ============================
%% save results to file
% ============================
% save uncertainties
filename = fullfile(stereo_results_directory, 'uncertainties.mat');
save(filename, 'individual_methods', 'unc_st_u', 'unc_st_v', 'unc_st_w', ...
    'Un_alpha1', 'Un_alpha2', 'Un_beta1', 'Un_beta2', 'tanalpha1', 'tanalpha2', 'tanbeta1', 'tanbeta2');

% save errors
filename = fullfile(stereo_results_directory, 'errors.mat');
save(filename, 'xgrid', 'ygrid', 'zgrid', 'Us', 'Vs', 'Ws', 'Ut', 'Vt', 'Wt', ...
            'err_u', 'err_v', 'err_w');

return;
% ============================
%% Plot RMS Uncertainty and Error fields (Spatial RMS Profile)
% ============================
fprintf('plotting RMS Uncertainty and Error fields (Spatial RMS Profile)\n');
% Good measurement cut off
ulx = 0.5;
uly = 0.5;
ulz = 1.5;

% initialize
err_u_rms = zeros(Sx, Sy);
err_v_rms = zeros(Sx, Sy);
err_w_rms = zeros(Sx, Sy);
unc_st_u_rms = {zeros(Sx, Sy); zeros(Sx, Sy); zeros(Sx, Sy)};
unc_st_v_rms = {zeros(Sx, Sy); zeros(Sx, Sy); zeros(Sx, Sy)};
unc_st_w_rms = {zeros(Sx, Sy); zeros(Sx, Sy); zeros(Sx, Sy)};

% loop through methods
for method_index = 1:num_individual_methods
    for i = 1:Sx
        for j = 1:Sy                
            % ============================
            % remove invalid measurements
            % ============================
            rms2derrx = squeeze(abs(err_u(i,j,:)));
            rms2dunx = squeeze(abs(unc_st_u{method_index}(i,j,:)));
            ix1 = find(rms2derrx>ulx);
            rms2derrx(ix1) = [];
            rms2dunx(ix1) = [];
            % Nx(i,j)=length(ix1);
            
            rms2derry = squeeze(abs(err_v(i,j,:)));
            rms2duny = squeeze(abs(unc_st_v{method_index}(i,j,:)));
            iy1 = find(rms2derry>uly);
            rms2derry(iy1) = [];
            rms2duny(iy1) = [];
            % Ny(i,j)=length(iy1);
            
            rms2derrz = squeeze(abs(err_w(i,j,:)));
            rms2dunz = squeeze(abs(unc_st_w{method_index}(i,j,:)));
            iz1 = find(rms2derrz>ulz);
            rms2derrz(iz1) = [];
            rms2dunz(iz1) = [];
            % Nz(i,j)=length(iz1);
        
            % ============================
            % calculate ensemble rms
            % ============================
            err_u_rms(i, j) = rms(rms2derrx);
            unc_st_u_rms{method_index}(i, j) = rms(rms2dunx);
            err_v_rms(i, j) = rms(rms2derry);
            unc_st_v_rms{method_index}(i, j) = rms(rms2duny);
            err_w_rms(i, j) = rms(rms2derrz);
            unc_st_w_rms{method_index}(i, j) = rms(rms2dunz);        
        end
    end

    % ============================
    % make plots
    % ============================
    figure;
    subplot(2,3,1); imagesc(err_u_rms); caxis([0 0.2]); title('Error, U')
    subplot(2,3,2); imagesc(err_v_rms); caxis([0 0.2]); title('Error, V')
    subplot(2,3,3); imagesc(err_w_rms); caxis([0 0.45]); title('Error, W')
    subplot(2,3,4); imagesc(unc_st_u_rms{method_index}); caxis([0 0.2]); title([individual_methods{method_index} ', U'])
    subplot(2,3,5); imagesc(unc_st_v_rms{method_index}); caxis([0 0.2]); title([individual_methods{method_index} ', V'])
    subplot(2,3,6); imagesc(unc_st_w_rms{method_index}); caxis([0 0.45]); title([individual_methods{method_index} ', W'])

    % save figures
    if save_figures
        save_figure_to_png_eps_fig(figure_save_directory, ['rms-spatial-' individual_methods{method_index}], [1, 0, 0]);
    end
end

% ============================
%% Plot RMS Uncertainty and Error fields (Temporal RMS Profile)
% ============================
% initialize
err_u_rms_fr = nans(1, num_frames);
err_v_rms_fr = nans(1, num_frames);
err_w_rms_fr = nans(1, num_frames);
unc_st_u_rms_fr = {nans(1, num_frames); nans(1, num_frames); nans(1, num_frames)};
unc_st_v_rms_fr = {nans(1, num_frames); nans(1, num_frames); nans(1, num_frames)};
unc_st_w_rms_fr = {nans(1, num_frames); nans(1, num_frames); nans(1, num_frames)};

for method_index = 1:num_individual_methods
    for frame_index=1:num_frames    
        % ============================
        % remove invalid measurements
        % ============================
        rex = squeeze(abs(err_u(:,:,frame_index)));
        rIMx = squeeze(abs(unc_st_u{method_index}(:,:,frame_index)));
        rex = rex(:);
        rIMx = rIMx(:);
        ix1 = find(rex>ulx);
        rex(ix1) = [];
        rIMx(ix1) = [];
        % Nx(i,j)=length(ix1);
        
        rey = squeeze(abs(err_v(:,:,frame_index)));
        rIMy = squeeze(abs(unc_st_v{method_index}(:,:,frame_index)));
        rey = rey(:);
        rIMy = rIMy(:);
        iy1 = find(rey>uly);
        rey(iy1) = [];
        rIMy(iy1) = [];
        % Ny(i,j)=length(iy1);
        
        rez = squeeze(abs(err_w(:,:,frame_index)));
        rIMz = squeeze(abs(unc_st_w{method_index}(:,:,frame_index)));
        rez = rez(:);
        rIMz = rIMz(:);
        iz1 = find(rez>ulz);
        rez(iz1) = [];
        rIMz(iz1) = [];
        % Nz(i,j)=length(iz1);
            
        % ============================
        % calculate rms
        % ============================
        err_u_rms_fr(frame_index) = rms(rex);
        unc_st_u_rms_fr{method_index}(frame_index) = rms(rIMx);
        
        err_v_rms_fr(frame_index) = rms(rey);
        unc_st_v_rms_fr{method_index}(frame_index) = rms(rIMy);
        
        err_w_rms_fr(frame_index) = rms(rez);
        unc_st_w_rms_fr{method_index}(frame_index) = rms(rIMz);
    end

    % ============================
    % make plots
    % ============================
    figure;hold on;
    plot(1:num_frames, err_u_rms_fr, 'r-o', 1:num_frames, unc_st_u_rms_fr{method_index}, 'r-*');
    plot(1:num_frames, err_v_rms_fr, 'b-o', 1:num_frames, unc_st_v_rms_fr{method_index}, 'b-*');
    plot(1:num_frames, err_w_rms_fr, 'g-o', 1:num_frames, unc_st_w_rms_fr{method_index}, 'g-*');
    hold off;
    title(individual_methods{method_index}); xlabel('Frame No.');
    ylabel('RMS error, uncertainty(pix/frame)');
    hleg=legend({'\sigma^{e}_{u}','\sigma^{CS}_{u}','\sigma^{e}_{v}','\sigma^{CS}_{v}','\sigma^{e}_{w}','\sigma^{CS}_{w}'}, ...
                'location', 'eastoutside');
    % set(hleg,'location','eastoutside');
    % set(gca,'FontSize',20);
    axis([0 50 0 0.5]); axis square;

    % save figures
    if save_figures
        save_figure_to_png_eps_fig(figure_save_directory, ['rms-temporal-' individual_methods{method_index}], [1, 0, 0]);
    end    
end

% ============================
% plot error and uncertainty pdfs
% ============================

% number of bins
Nc=100;
% max uncertainty
ul=0.5;

% bins
bins_u1 = linspace(0, ul, Nc);
bins_v1 = linspace(0, ul, Nc);
bins_u2 = linspace(0, ul, Nc);
bins_v2 = linspace(0, ul, Nc);

% initialize
N_u1 = nans(num_individual_methods, Nc);
N_v1 = nans(num_individual_methods, Nc);
N_u2 = nans(num_individual_methods, Nc);
N_v2 = nans(num_individual_methods, Nc);

% loop through methods
for method_index = 1:num_individual_methods
    % calculate pdfs
    N_u1(method_index, :) = histc(unc1_u{method_index}(:), bins_u1);
    N_v1(method_index, :) = histc(unc1_v{method_index}(:), bins_v1);
    N_u2(method_index, :) = histc(unc2_u{method_index}(:), bins_u2);
    N_v2(method_index, :) = histc(unc2_v{method_index}(:), bins_v2);

    % plot pdfs
    figure;    
    plot(bins_u1, N_u1(method_index, :), 'r-o', bins_v1, N_v1(method_index, :), 'b-o', ...
            bins_u2, N_u2(method_index, :), 'r-+', bins_v2, N_v2(method_index, :), 'b-+');
    title('Planar Correlation Uncertainties');
    legend('U1','V1','U2','V2');
    if save_figures
        % print(gcf,'-dpng',fullfile(figure_save_directory,'planarhist.png'),'-r300');
        save_figure_to_png_eps_fig(figure_save_directory, ['histogram-planar-' individual_methods{method_index}], [1, 0, 0]);
    end
end

% ============================
%% Calculate coverage
% ============================
% initialize
cov_u = nans(1, num_individual_methods);
cov_v = nans(1, num_individual_methods);
cov_w = nans(1, num_individual_methods);

% number of bins
Nb = 30;

% bins
bins_err_u = linspace(-ulx, ulx, 2*Nb-1);
bins_err_v = linspace(-uly, uly, 2*Nb-1);
bins_err_w = linspace(-ulz, ulz, 2*Nb-1);
bins_unc_st_u = linspace(0, ulx, Nb);
bins_unc_st_v = linspace(0, uly, Nb);
bins_unc_st_w = linspace(0, ulz, Nb);

% initialize
N_err_u = nans(1, 2*Nb-1);
N_err_v = nans(1, 2*Nb-1);
N_err_w = nans(1, 2*Nb-1);
N_unc_st_u = nans(num_individual_methods, Nb);
N_unc_st_v = nans(num_individual_methods, Nb);
N_unc_st_w = nans(num_individual_methods, Nb);

% loop through methods
for method_index = 1:num_individual_methods
    % ============================
    % aggregate measurements
    % ============================
    err_u_temp = err_u(:);
    err_v_temp = err_v(:);
    err_w_temp = err_w(:);

    unc_st_u_temp = unc_st_u{method_index}(:);
    unc_st_v_temp = unc_st_v{method_index}(:);
    unc_st_w_temp = unc_st_w{method_index}(:);

    % ============================
    % eliminate invalid measurements
    % ============================
    inx = find(abs(err_u_temp(:)) > ulx);
    iny = find(abs(err_v_temp(:)) > uly);
    inz = find(abs(err_w_temp(:)) > ulz);

    err_u_temp(inx) = [];
    err_v_temp(iny) = [];
    err_w_temp(inz) = [];
    unc_st_u_temp(inx) = [];
    unc_st_v_temp(iny) = [];
    unc_st_w_temp(inz) = [];

    % ============================
    % calculate coverage
    % ============================
    cnt1=0;cnt2=0;cnt3=0;
    l1=length(err_u_temp);
    l2=length(err_v_temp);
    l3=length(err_w_temp);
    for m1=1:l1
        if abs(err_u_temp(m1))<=unc_st_u_temp(m1)
            cnt1=cnt1+1;
        end    
    end

    for m2=1:l2
        if abs(err_v_temp(m2))<=unc_st_v_temp(m2)
            cnt2=cnt2+1;
        end    
    end

    for m3=1:l3    
        if abs(err_w_temp(m3))<=unc_st_w_temp(m3)
            cnt3=cnt3+1;
        end    
    end

    cov_u(method_index) = 100 * cnt1/l1;
    cov_v(method_index) = 100 * cnt2/l2;
    cov_w(method_index) = 100 * cnt3/l3;
    fprintf('\n Coverage with angle uncertainty \n Ucov=%2.2f Vcov=%2.2f Wcov=%2.2f \n', ...
            [cov_u(method_index) cov_v(method_index) cov_w(method_index)]);

    set(0,'DefaultAxesFontName', 'Times New Roman');

    % ============================
    % calculate pdfs
    % ============================
    N_err_u = histc(err_u_temp, bins_err_u) / length(err_u_temp);
    N_unc_st_u(method_index, :) = histc(unc_st_u_temp, bins_unc_st_u) / length(err_u_temp);

    N_err_v = histc(err_v_temp, bins_err_v) / length(err_v_temp);
    N_unc_st_v(method_index, :) = histc(unc_st_v_temp, bins_unc_st_v) / length(err_v_temp);

    N_err_w = histc(err_w_temp, bins_err_w) / length(err_w_temp);
    N_unc_st_w(method_index, :) = histc(unc_st_w_temp, bins_unc_st_w) / length(err_w_temp);

    % ============================
    % Plot pdfs
    % ============================

    figure;
    
    % ============================
    % U uncertainty Histogram
    % ============================
    subplot(1, 3, 1)
    hold on;
    plot(bins_err_u, N_err_u, '-k');
    plot(bins_unc_st_u, N_unc_st_u(method_index, :), 'r-');
    % HL1=fill(bins_err_u,N_err_u,'-k');
    % HL2=fill(bins_unc_st_u,N_unc_st_u,'r-');
    % set(HL1,'facealpha',.5);
    % set(HL2,'facealpha',.5);
    yl = ylim(gca);
    ll = linspace(0, max(yl(2)), 10);
    sigu = rms(err_u_temp(abs(err_u_temp) <= ulx)) .* ones(length(ll),1);
    mixx1 = rms(unc_st_u_temp(unc_st_u_temp <= ulx)) .* ones(length(ll),1);
    plot(sigu, ll, 'k--', mixx1, ll, 'r--'); 

    grid on
    box on 
    set(gcf,'color','white');
    legend('e_{u}', ['\sigma_{u}^{',individual_methods{method_index},'}']);
    axis([-ulx ulx 0 inf]);
    % title(['U Uncertainty Histogram coverage= ',num2str(cov_u,'%02.02f')]);
    title('U');
    xlabel('Error/Uncertainty (pix.)'); ylabel('No. of vectors');

    % ============================
    % V uncertainty Histogram
    % ============================
    subplot(1, 3, 2)
    hold on;
    plot(bins_err_v, N_err_v, '-k');
    plot(bins_unc_st_v, N_unc_st_v(method_index, :), 'b-');
    yl = ylim(gca);
    ll = linspace(0, max(yl(2)), 10);
    sigv = rms(err_v_temp(abs(err_v_temp) <= uly)) .* ones(length(ll),1);
    miyy1 = rms(unc_st_v_temp(unc_st_v_temp <= uly)) .* ones(length(ll),1);
    plot(sigv, ll, 'k--', miyy1, ll, 'b--'); 

    grid on
    box on 
    set(gcf, 'color', 'white');
    legend('e_{v}', ['\sigma_{v}^{',individual_methods{method_index},'}']);
    axis([-uly uly 0 inf]);
    % title(['V uncertainty histogram coverage= ',num2str(cov_v,'%02.02f')]);
    title('V');
    xlabel('Error/Uncertainty (pix.)'); ylabel('No. of vectors');

    % ============================
    % W uncertainty Histogram
    % ============================
    subplot(1, 3, 3)
    hold on;
    plot(bins_err_w, N_err_w, '-k');
    plot(bins_unc_st_w, N_unc_st_w(method_index, :), 'g-');
    yl = ylim(gca);
    ll = linspace(0, max(yl(2)), 10);
    sigw = rms(err_w_temp(abs(err_w_temp) <= ulz)) .* ones(length(ll),1);
    mizz1 = rms(unc_st_w_temp(unc_st_w_temp <= ulz)) .* ones(length(ll),1);
    plot(sigw, ll, 'k--', mizz1, ll, 'g--'); 

    grid on
    box on
    set(gcf, 'color', 'white');
    legend('e_{w}', ['\sigma_{w}^{',individual_methods{method_index},'}']);
    axis([-ulz ulz 0 inf]);
    % title(['W uncertainty histogram coverage= ',num2str(cov_w,'%02.02f')]);
    title('W');
    xlabel('Error/Uncertainty (pix.)'); ylabel('No. of vectors');

    set(gcf, 'resize', 'off');
    set(gcf, 'Position', [80         364        1002         306]);
    if save_figures
        save_figure_to_png_eps_fig(figure_save_directory, ['pdf-stereo-' individual_methods{method_index}], [1, 0, 0])
        % export_fig(gcf,fullfile([figure_save_directory,'Whist.png']),'-painters','-r360');
    end

    fprintf('\n RMS error & Uncertainty for x,y,z');
    fprintf('\n ex,ux,ey,uy,ez,uz \n %0.3f %0.3f %0.3f %0.3f %0.3f %0.3f \n', [sigu(1) mixx1(1) sigv(1) miyy1(1) sigw(1) mizz1(1)]);
end

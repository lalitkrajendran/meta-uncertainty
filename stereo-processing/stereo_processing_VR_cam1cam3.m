% Script to load and process 2d images, and perform stereo 
% reconstruction for vortex ring cam1 cam3

clc;
clear;
close all;

restoredefaultpath;
addpath(genpath('/scratch/shannon/c/aether/Projects/BOS/general-codes/matlab-codes/'));
addpath('../prana-uncertainty-average-dc-new-im-2-save-deform-cs/');
% addpath /home/shannon/a/bhattac3/Stereo_Uncertainty_Work/Codes/prana_new_version_10_28_2015/prana/;
setup_default_settings;

% ====================================
%% read/write settings
% ====================================
% top level read directory
top_read_directory = '/scratch/shannon/a/aether/Projects/PIV_PTV_Uncertainty_Quantification/Stereo_uncertainty_work/';
% top level write directory
top_write_directory = '/scratch/shannon/a/aether/Projects/PIV_PTV_Uncertainty_Quantification/uncertainty-combination/analysis/results/stereo-dataset/cam13/With_selfcal/';
mkdir_c(top_write_directory);

% Load Calibration job file;
% caljobfile=load('/home/shannon/a/bhattac3/Stereo_Uncertainty_Work/New_Stereo_uncertainty_tests/Vortex_Ring/Prana_process/jobs/cam13_afterselfcalibration_job2.mat');
caljobfile=load(fullfile(top_read_directory, 'Jobs', 'cam13_afterselfcalibration_job2.mat'));

% Load 2d job file 
% pranajob=load('/home/shannon/a/bhattac3/Stereo_Uncertainty_Work/New_Stereo_uncertainty_tests/Vortex_Ring/Prana_process/results/cam13/With_selfcal/vortexring13_2djob.mat');
pranajob=load(fullfile(top_read_directory, 'Results', 'Prana_process/results/cam13/With_selfcal/vortexring13_2djob.mat'));

% ===================================
% modify job file
% ===================================
rectype='Willert';
planarjob=pranajob.Data;
caldata=caljobfile.datasave.caljob;

% change read directory
% planarjob.imdirec = '/scratch/shannon/a/aether/Projects/PIV_PTV_Uncertainty_Quantification/Stereo_uncertainty_work/Images/Vortex_Ring_Data/E_Particle_Images/Camera_01/';
planarjob.imdirec = fullfile(top_read_directory, 'Images', 'Vortex_Ring_Data/E_Particle_Images/Camera_01/');
planarjob.imdirec2 = fullfile(top_read_directory, 'Images', 'Vortex_Ring_Data/E_Particle_Images/Camera_03/');

% change write directory
% planarjob.outdirec = '/scratch/shannon/a/aether/Projects/PIV_PTV_Uncertainty_Quantification/Stereo_uncertainty_work/Results/Prana_process/results/cam13/With_selfcal/';
planarjob.outdirec = top_write_directory;

pulsesep = str2double(planarjob.wrsep)*(1e-6);
planarjob.wrsep = '1';
planarjob.wrmag = '1';

% temporary change to quickly troubleshoot results
% planarjob.imfend = planarjob.imfstart;
% planarjob.par = '0';

% turn on saving deformed images
planarjob.SaveIMdeform = '1';

% =================================================================
% turn on uncertainty calculation for the last and penultimate pass
% =================================================================
% set uncertainty flags for individual passes
for pass_index=0:str2double(planarjob.passes)
    if pass_index < str2double(planarjob.passes)
        uncertaintyestimate = '0';
    else
        uncertaintyestimate = '1';
    end
    planarjob.(['PIV' num2str(pass_index)]).write = '1';
    planarjob.(['PIV' num2str(pass_index)]).uncertaintyestimate = uncertaintyestimate;
    planarjob.(['PIV' num2str(pass_index)]).ppruncertainty = uncertaintyestimate;
    planarjob.(['PIV' num2str(pass_index)]).miuncertainty = uncertaintyestimate;
    planarjob.(['PIV' num2str(pass_index)]).imuncertainty = uncertaintyestimate;
    planarjob.(['PIV' num2str(pass_index)]).mcuncertainty = uncertaintyestimate;
    planarjob.(['PIV' num2str(pass_index)]).csuncertainty = uncertaintyestimate;
end

% turn on validation for last pass
planarjob.PIV4.val = '1';
planarjob.PIV4.uod = '1';
planarjob.PIV4.thresh = '1';
planarjob.PIV4.velsmooth = '1';

% ============================
%% Dewarp the camera images
% ============================
fprintf('Processing for Geometric Reconstruction... \n');
% [dewarpdirlist,dewarp_grid,wil_scaling]=imagedewarp(caldata,'Willert',planarjob);
% [Xg,Yg] = meshgrid(-34:.45:20,-27:.45:27);
[Xg,Yg] = meshgrid(linspace(-34,20,1024), linspace(-27,27,1024));
[dewarpdirlist,dewarp_grid,wil_scaling]=imagedewarp_predefined_grid(caldata,'Willert',planarjob,[],Xg,Yg,8);
diroutlist.dewarpdirlist=dewarpdirlist;

% ============================
%% 2D processing
% ============================
for camera_index = 1:2
    % extract job file
    job = planarjob;
    job.imdirec = dewarpdirlist.(['dewarpdir' num2str(camera_index)]);
    job.imdirec2 = jobfile.imdirec;
end

% ============================
%% 2D processing for camera1
% ============================
job1=planarjob;
job1.imdirec=dewarpdirlist.dewarpdir1;
job1.imbase=planarjob.imbase;
job1.imdirec2=dewarpdirlist.dewarpdir1;
job1.imbase2=planarjob.imbase;
% By creating a variable "cam1dir" I know only have to update one
% location if I want to change the form of this directory.  The
% previous way it was written I would have to update multiple locations
% increasing the chance I would miss one of them.
cam1dir = fullfile(job1.outdirec,rectype,['Camera',num2str(caldata.camnumber(1)),filesep], 'vectors');
mkdir_c(cam1dir);

job1.outdirec=cam1dir;
diroutlist.willert2dcam1=job1.outdirec;
fprintf(['\nProcessing Planar Fields for Camera:',num2str(caldata.camnumber(1)),'\n']);

pranaPIVcode(job1);

% ============================
%% 2D processing for camera2
% ============================
job2=planarjob;
job2.imdirec=dewarpdirlist.dewarpdir2;
job2.imbase=planarjob.imbase2;
job2.imdirec2=dewarpdirlist.dewarpdir2;
job2.imbase2=planarjob.imbase2;

cam2dir = fullfile(job2.outdirec,rectype,['Camera',num2str(caldata.camnumber(2)),filesep], 'vectors');
mkdir_c(cam2dir);

job2.outdirec=cam2dir;
diroutlist.willert2dcam2=job2.outdirec;
fprintf(['\nProcessing Planar Fields for Camera:',num2str(caldata.camnumber(2)),'\n']);
pranaPIVcode(job2);

% ============================
%% stereo reconstruction
% ============================
stereodir = fullfile(planarjob.outdirec, rectype, ['Camera',num2str(caldata.camnumber(1)),'Camera',num2str(caldata.camnumber(2)),'_3Cfields',filesep], 'vectors');
mkdir_c(stereodir);
diroutlist.willert3cfields = stereodir;

fprintf('Doing Geometric Stereo Reconstructions.... \n')
willert_vec_reconstruct(diroutlist,caldata,dewarp_grid,wil_scaling,pulsesep);
scaling.wil = wil_scaling;

%%  save jobfile
save(fullfile(top_write_directory, 'cam13_VR_prana_fulljob_withselfcal.mat'), 'caljobfile', 'planarjob', 'diroutlist', 'scaling', 'dewarp_grid', 'job1', 'job2');



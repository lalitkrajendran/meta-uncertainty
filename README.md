This folder contains the codes used for the uncertainty combination/meta-model methodology outlined in the following paper.

# Paper
Rajendran, L., Bhattacharya, S., Bane, S., & Vlachos, P. (2021). Meta-Uncertainty for Particle Image Velocimetry. Measurement Science and Technology. https://doi.org/10.1088/1361-6501/abf44f

# Code Organization
Codes are in the src folder and organized into the following subfolders. 
Many codes call PRANA. Clone from: https://github.itap.purdue.edu/lrajendr/prana and add to path.

## meta-uncertainty-code-package
codes to perform a sample meta-uncertainty calculation

- sample_script.m: script implements all steps of the meta uncertainty calculation for a single interrogation window

required data also provided in this same folder
- sample-displacements.mat: displacements from prana processing
- jobfile.mat: processing details
- im1-def.mat: deformed image 1
- im2-def.mat: deformed image 2

Code saves to sample-result.mat 

## planar-processing
codes to process the planar piv images, and post-process the results to calculate the 

* batch_prana_process_images.m: to process the images
* batch_display_vector_fields.m: to display the vectors from processing
* batch_calculate_physical_properties.m: to calculate global properties like velocity gradients for assimilating data
* batch_calculate_errors_02.m: to calculate errors


* batch_run_resampling_monte_carlo_individual_02.m: to run the resampling
* batch_calculate_weights_individual_03.m: to calculate weights
* batch_calculate_error_uncertainty_statistics_individual_02.m: to calculate the error and uncertainty statistics

* plot_weighting_based_on_ratio_change_04.m: plot results in assimilated form
* plot_weighting_based_on_ratio_change_panel.m: plot results in individual panel form
* plot_effect_of_resampling.m: violin plots to effect of resampling methods on uncertainty and snr

## stereo-processing
set of codes to process stereo piv images and vectors

- stereo_processing_VR_cam1cam3.m: process stereo images
- calculate_stereo_uncertainty_vortex_ring.m: calculate the stereo uncertainty using Sayantan's codes in stereo_uncertainty_codes_packaged

- batch_run_resampling_stereo_monte_carlo_02.m: codes to run resampling 
- batch_calculate_weights_stereo_3c_02: calculate weights from resampling results
- calculate_stereo_error_uncertainty_statistics_monte_carlo_05.m: calculate rms error/unc and pdfs to make the plots

- plot_stereo_error_uncertainty_statistics_monte_carlo_08.m: make the plots in the paper

## general-codes
set of functions for processing images and plotting results that are common to both the planar and stereo processing

all function files contain documentation explaining their operation

# Results

## planar-dataset
results for the planar processing

- experiment-new: contains vector fields from planar piv processing for each dataset
- monte-carlo: contains resampling results for all datasets as well as final results and figures used in the paper

## stereo-dataset
results for stereo processing. each folder contains both the processed vector fields and the resampling results for the respective configurations

- cam13: processing for cam1, cam3 configuration. Refer Sayantan's stereo uncertainty paper for experimental arrangement.
- cam24: processing for cam2, cam4 configuration. Refer Sayantan's stereo uncertainty paper for experimental arrangement.







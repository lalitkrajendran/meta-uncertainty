clear 
close all
clc

%%
read_directory = '/scratch/shannon/a/aether/Projects/PIV_PTV_Uncertainty_Quantification/Planar_uncertainty_work/Images/experiment/Jetdata/';

[files, num_files] = get_directory_listing(read_directory, 'B*.tif');

write_directory = '/scratch/shannon/a/aether/Projects/PIV_PTV_Uncertainty_Quantification/Planar_uncertainty_work/Images/experiment/Jetdata_mod/';
mkdir_c(write_directory);

%%
for file_index = 1:num_files
    fprintf('file_index: %d\n', file_index);
%     new_file_name = [files(file_index).name(1:end-3), 'tif'];
    new_file_name = sprintf('im_%04d.tif', file_index);
    copyfile(fullfile(read_directory, files(file_index).name), fullfile(write_directory, new_file_name));
end

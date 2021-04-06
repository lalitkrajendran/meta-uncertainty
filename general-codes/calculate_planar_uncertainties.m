% function [unc_planar, snr_metric, uncertainty2D] = calculate_planar_uncertainties(im1, im2, pass_settings, uncertainty_methods)
function [unc_planar, snr_metric, U, V, uncertainty2D] = calculate_planar_uncertainties(im1, im2, pass_settings, uncertainty_methods)

    % number of uncertainty methods
    num_methods = numel(uncertainty_methods);

    % calculate image size
    [image_height, image_width] = size(im1);

    % extract window size
    [window_size_x, window_size_y] = extract_window_size(pass_settings);
    window_size = [window_size_x, window_size_y];

    % extract window resolution
    [window_resolution_x, window_resolution_y] = extract_window_resolution(pass_settings);
    window_resolution = [window_resolution_x, window_resolution_y];

    % extract correlation diameter
    str = strsplit(pass_settings.RPCd, ',');
    D = [str2double(str{1}) str2double(str{2})];

    % extract uncertainty flags
    uncertainty_flags = extract_uncertainty_flags(pass_settings);

    % create cell to hold uncertainties
    unc_planar = cell(1, num_methods);
    
    % loop through all uncertainty methods
    for method_index = 1:num_methods
        % create structure to hold uncertainties for the current method
        unc_planar{method_index} = struct;
        % assign name of the current method to the struct
        unc_planar{method_index}.name = uncertainty_methods{method_index};
            
        % calculate uncertainty using IM        
        if strcmp(uncertainty_methods{method_index}, 'IM')
            % [sigma_x, sigma_y, ~, ~] = original_particle_disparity(im1, im2, image_width/2, image_height/2, 0, 0, window_resolution(1), 0);
            % [sigma_x, sigma_y, ~, ~] = original_particle_disparity(im1, im2, image_width/2, image_height/2, 0, 0, window_resolution(1));
            [sigma_x, sigma_y, ~, ~] = original_particle_disparity(im1, im2, image_width/2, image_height/2, window_resolution(1));
        
        % calculate uncertainty using MC
        elseif strcmp(uncertainty_methods{method_index}, 'MC')            
            % [~, ~, U, V, ~, ~, ~, uncertainty2D, SNRmetric] = PIVwindowed(im1, im2, 'SCC', window_size, [window_resolution; window_resolution], 0, [2.8, 2.8], 1, 3, 0, 0, 0, image_width/2, image_height/2, uncertainty_flags, 1, 0, 0);            
            % [Xc,Yc,Uc,Vc,Cc,Dc,Cp,uncertainty2D,SNRmetric]=PIVwindowed(im1d,im2d,Corr{e},Wsize(e,:),Wres(:, :, e),0,D(e,:),Zeromean(e),Peaklocator(e),find_extrapeaks,frac_filt(e),saveplane(e),X(Eval>=0),Y(Eval>=0),uncertainty,e);
            [~, ~, U, V, ~, ~, ~, uncertainty2D, SNRmetric] = PIVwindowed(im1, im2, pass_settings.corr, window_size, [window_resolution; window_resolution], 0, D, ...
                                str2double(pass_settings.zeromean), str2double(pass_settings.peaklocator), str2double(pass_settings.savepeakinfo), ...
                                str2double(pass_settings.frac_filt), str2double(pass_settings.savepeakinfo), image_width/2, image_height/2, uncertainty_flags, 1, 0, 0);            

            sigma_x = sqrt(uncertainty2D.biasx.^2 + (uncertainty2D.Ixx.^2) ./ uncertainty2D.Neff);
            sigma_y = sqrt(uncertainty2D.biasy.^2 + (uncertainty2D.Iyy.^2) ./ uncertainty2D.Neff);
            % sigma_x = uncertainty2D.Ixx ./ sqrt(uncertainty2D.Neff);
            % sigma_y = uncertainty2D.Iyy ./ sqrt(uncertainty2D.Neff);
            snr_metric = extract_snr_metric(SNRmetric, 1, 1);
            
        % calculate uncertainty using CS
        elseif strcmp(uncertainty_methods{method_index}, 'CS')
            [sigma_x, sigma_y] = correlation_statistics(im1, im2, window_size, [window_resolution; window_resolution], window_size(1)/2, window_size(2)/2);        
        end

        % allocate uncertainties
        unc_planar{method_index}.x = sigma_x;
        unc_planar{method_index}.y = sigma_y;
    end
    
end

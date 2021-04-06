function arr_valid_all = extract_valid_resampled_elements(arr, valid_trials, num_resampling_trials)
% Function to extract valid resampled elements corresponding to each valid trial
%
% INPUTS:
% arr: array containing elements to extracted
% valid_trials: array of indices/trials that are considered valid
% num_resampling_trials: number of resampling trials for each global trial
% 
% OUTPUTS:
% arr_valid_all: array containing valid elements only
%
% AUTHOR:
% Lalit Rajendran (lrajendr@purdue.edu)
%
% DATE:
% 2020/06/19
    
    arr_valid_all = [];
    for trial_index = valid_trials
        % start index
        start_index = (trial_index - 1) * num_resampling_trials + 1;
        % stop index
        stop_index = trial_index * num_resampling_trials;
        % extract valid resampled values
        arr_valid_current = arr(start_index:stop_index);
        % aggregate
        arr_valid_all = [arr_valid_all, arr_valid_current];            
    end
end
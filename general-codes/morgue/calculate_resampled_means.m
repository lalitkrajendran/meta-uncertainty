function arr_mean = calculate_resampled_means(arr, num_resampling_trials)
% Function to calculate means for each set of resampling trials
%
% INPUTS:
% arr: array containing values whose mean is to be computed
% num_resampling_trials: number of resampling trials for each global trial
% 
% OUTPUTS:
% arr_mean: array of mean values that correspond to each set of resampling trials
%
% AUTHOR:
% Lalit Rajendran (lrajendr@purdue.edu)
%
% DATE:
% 2020/06/19
    
    num_trials = numel(arr)/num_resampling_trials;
    
    for trial_index = 1:num_trials
        % start index
        start_index = (trial_index - 1) * num_resampling_trials + 1;
        % stop index
        stop_index = trial_index * num_resampling_trials;
        % calculate mean
        arr_mean(trial_index) = mean(arr(start_index:stop_index), 'omitnan');                                
    end 
end
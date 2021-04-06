function uncertainty_flags = extract_uncertainty_flags(pass_settings)
% Function to extract uncertainty flags from the pass settings array
%
% INPUTS:
% pass_settings: jobfile settings for this pass
%
% OUTPUTS:
% uncertainty_flags: structure of integer 0/1 flags for all uncertainty schemes
%
% AUTHORS:
% Lalit Rajendran (lrajendr@purdue.edu)
%
% DATE:
% 06/12/2020

    % create structure
    uncertainty_flags = struct;    
    % extract field names
    field_names = fieldnames(pass_settings);
    % calculate number of field names
    num_fields = numel(field_names);

    % loop through fields
    for field_index = 1:num_fields
        % extract name of the current field
        current_field_name = field_names{field_index};
        % extract value of the current field
        current_value = pass_settings.(current_field_name);

        % if the field contains the word 'uncertainty', then extract its value
        if contains(current_field_name, 'uncertainty')
            % if the flag is a character, then convert to integer
            if ischar(current_value)
                uncertainty_flags.(current_field_name) = str2double(current_value);
            % else, just copy the value onto the structure
            else
                uncertainty_flags.(current_field_name) = current_value;
            end
        end
    end
    
end


function mi = calculate_mutual_information(x1, x2, bins1, bins2)
% Function to calculate mutual information between two vectors of
% observations. The mutual information is defined as the difference
% between the sum of the individual entropies and the joint entropy.
%
% INPUTS:
% x1, x2: vectors of observations
% bins1, bins2: bins for histogram calculations
%
% OUTPUTS:
% mi: mutual information between the two variables
%
% AUTHOR
% Lalit Rajendran (lrajendr@purdue.edu)
%
% DATE:
% 2020/06/08

    % calculate self entropy
    e1 = calculate_self_entropy(x1, bins1);
    e2 = calculate_self_entropy(x2, bins2);

    % calculate joint entropy
    e12 = calculate_joint_entropy(x1, x2, bins1, bins2);

    % calculate mutual information
    mi = e1 + e2 - e12;
    
end
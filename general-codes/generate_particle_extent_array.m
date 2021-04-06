function [r, c] = generate_particle_extent_array(xp, yp, dp)
% Function to generate array locations on an image corresponding to a particle
%
% INPUTS:
% xp, yp: particle centroids (pix.)
% dp: particle diameter (pix.)
%
% OUTPUTS:
% r, c: row and column locations corresponding to particle intensity
%
% AUTHOR:
% Lalit Rajendran (lrajendr@purdue.edu)

    % round value
    d = round(dp);

    % generate extent
    particle_extent = (-d:d);

    % calculate array locations on the image corresponding to the particle 
    r = particle_extent + round(yp);
    c = particle_extent + round(xp);

end
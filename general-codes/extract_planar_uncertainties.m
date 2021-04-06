function unc_planar = extract_planar_uncertainties(uncertainty2D, r, c)
% Function to extract planar uncertainties for a 
% given point in the FOV
%
% INPUTS:
% uncertainty2D: structure from prana results containing uncertainties
% r, c: row and column indices of the grid point
%
% OUTPUTS:
% unc_planar: structure containing planar uncertainties for this grid point
%
% AUTHORS:
% Lalit Rajendran (lrajendr@purdue.edu)
%
% DATE:
% 2020/04/28

    unc_planar = struct;
    % image matching
    unc_planar.imx = uncertainty2D.Uimx(r, c);
    unc_planar.imy = uncertainty2D.Uimy(r, c);
    % moment of correlation
    unc_planar.mcx = uncertainty2D.MCx(r, c);
    unc_planar.mcy = uncertainty2D.MCy(r, c);
    % correlation statistics
    unc_planar.csx = uncertainty2D.Ucsx(r, c);
    unc_planar.csy = uncertainty2D.Ucsy(r, c);
end
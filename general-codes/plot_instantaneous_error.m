function plot_instantaneous_error(X, Y, U, V, Ut, Vt, err_U, err_V)
% Function to plot contours of instantaneous errors along U and V
%
% INPUTS:
% X, Y: co-ordinate grid
% U, V: displacements along X and Y
% Ut, Vt: true, interpolated solution along X and Y on the same grid points
% err_U, err_V: errors along X and Y
%
% OUTPUTS:
% None
%
% AUTHORS:
% Lalit Rajendran (lrajendr@purdue.edu)
%
% DATE:
% 2020/06/19

    figure
    
    % processing
    subplot(2, 3, 1)
    contourf(X, Y, U)
    axis equal
    colorbar
    % get colorbar limits
    clim = caxis;
    title('U, Processed')
    
    % true solution, interpolated
    subplot(2, 3, 2)
    contourf(X, Y, Ut)
    caxis(clim)
    axis equal
    title('U, True, Interp')
    
    % error
    subplot(2, 3, 3)
    % contourf(X, Y, abs(err_U(:, :, 1))), caxis(clim)
    contourf(X, Y, abs(err_U))
    caxis(clim)
    axis equal
    title('U, Error')

    % prana
    subplot(2, 3, 4)
    contourf(X, Y, V)
    axis equal
    colorbar
    % get colorbar limits
    clim = caxis;
    title('V, Processed')
    
    % true solution, interpolated
    subplot(2, 3, 5)
    contourf(X, Y, Vt)
    caxis(clim)
    axis equal
    title('V, True, Interp')
    
    % error
    subplot(2, 3, 6)
    % contourf(X, Y, abs(err_U(:, :, 1))), caxis(clim)
    contourf(X, Y, abs(err_V))
    caxis(clim)
    axis equal
    title('V, Error')
end
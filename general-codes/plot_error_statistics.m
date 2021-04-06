function plot_error_statistics(X, Y, bias_U, bias_V, random_U, random_V, total_U, total_V)
% Function to plot contours of error statistics along U and V
%
% INPUTS:
% X, Y: co-ordinate grid
% bias_U, bias_V: bias error along X and Y
% random_U, random_V: random error along X and Y
% total_U, total_V: total error along X and Y
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
    
    % bias, U
    subplot(2, 3, 1)
    contourf(X, Y, abs(bias_U))
    axis equal
    colorbar
    % get colorbar limits
    clim = caxis;
    title('U, bias')
    
    % random, U
    subplot(2, 3, 2)
    contourf(X, Y, random_U)
    caxis(clim)
    axis equal
    title('U, Random')
    
    % total, U
    subplot(2, 3, 3)
    contourf(X, Y, total_U)
    caxis(clim)
    axis equal
    title('U, Total')

    % bias, V
    subplot(2, 3, 4)
    contourf(X, Y, abs(bias_V))
    axis equal
    colorbar
    % get colorbar limits
    clim = caxis;
    title('V, bias')
    
    % random, V
    subplot(2, 3, 5)
    contourf(X, Y, random_V)
    caxis(clim)
    axis equal
    title('V, Random')
    
    % total, V
    subplot(2, 3, 6)
    contourf(X, Y, total_V)
    caxis(clim)
    axis equal
    title('V, Total')
end
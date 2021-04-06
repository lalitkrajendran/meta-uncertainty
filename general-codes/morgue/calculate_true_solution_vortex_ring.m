function [true_solution, c, r] = calculate_true_solution_vortex_ring(true_solution, X, Y)
% This function converts the true solution for the vortex ring datasets for
% error analysis with Prana. It convert displacements to pixel/frame and
% interpolates it onto the processing grid.

    %% define parameters
    N=50;
    fstart=1;
    foffsetprana=24;
    foffsetdavis=1;

    spf=1;
    xscale=1/11.9690;%stereojob.datasave.scaling.wil.xscale;
    yscale=1/11.9653;%stereojob.datasave.scaling.wil.yscale;
    xorigin=576.8623;
    yorigin=395.8564;

    % true solution grid in object space
%     [Xt,Yt] = meshgrid(-34:0.45:20,-27:0.45:27);
    
    %% map true solution grid to pixel grid using magnification and
    % origin shift
    true_solution.Xt_scaled = true_solution.Xt./xscale+xorigin;
    true_solution.Yt_scaled = true_solution.Yt./yscale+yorigin;
    
    %% get limits of true solution grid in pixel coordinates
    true_solution.xmin = min(true_solution.Xt_scaled(:));
    true_solution.xmax = max(true_solution.Xt_scaled(:));
    true_solution.ymin = min(true_solution.Yt_scaled(:));
    true_solution.ymax = max(true_solution.Yt_scaled(:));
    
    %% calculate scaled velocities    
    true_solution.Ut_scaled = (spf/xscale) .* true_solution.Ut;
    true_solution.Vt_scaled = (spf/xscale) .* true_solution.Vt;
    
    %% find vector locations in the processing results for which
    % true solution is available
    
    % extract co-ordinate arrays for the processing grid
    xvec = X(1, :);
    yvec = Y(:, 1);
    
    % find vector locations in the processing results for which
    % true solution is available
    c = find(xvec > true_solution.xmin & xvec < true_solution.xmax);
    r = find(yvec > true_solution.ymin & yvec < true_solution.ymax);

end
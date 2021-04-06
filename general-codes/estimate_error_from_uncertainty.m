function err_est = estimate_error_from_uncertainty(model, avg, unc)
% Function to estimate error from uncertainty, for a gaussian or a log-normal distribution
%
% INPUTS:
% unc: uncertainty
% model: model for the pdf 'gaussian' vs 'lognormal'
%
% OUTPUTS:
% err_est: estimated error
%
% AUTHOR:
% Lalit Rajendran (lrajendr@purdue.edu)

    % number of elements
    N = numel(unc);
    % number of trials for each grid point
    num_trials = 1e3;
    % replicate arrays to generate 1000 random numbers for each uncertainty vector
    avg = repmat(avg, 1, num_trials);
    unc = repmat(unc, 1, num_trials);

    rng();
    % estimate error from a Gaussian pdf
    if strcmp(model, 'gaussian')
        % err_est = randn(N, num_trials) .* unc;
        err_est = randn(N, num_trials) .* unc + avg;
    
    % estimate error from a lognormal pdf
    elseif strcmp(model, 'lognormal')
        % estimate sigma from uncertainty
        % unc
        % y = 0.5 + sqrt(0.25 + unc.^2)
        % sig = sqrt(log(y))
        % err_est = lognrnd(0, unc, N, 1) - 1;

        % construct variables for lognormal distribution
        mu = log(avg.^2 ./ sqrt(avg.^2 + unc.^2));
        sig = sqrt(log(1 + unc.^2 ./ avg.^2));
        err_est = lognrnd(mu, sig);
    
    % estimate error from an exponential pdf
    elseif strcmp(model, 'exp')
        err_est = exprnd(unc);
    
    % estimate error from a rayleigh pdf
    elseif strcmp(model, 'rayleigh')
        % calculate scale parameter
        b = sqrt(2 / (4 - pi)) * unc/sqrt(2);
        
        % estimate error
        err_est = raylrnd(b);
    end       
    
end
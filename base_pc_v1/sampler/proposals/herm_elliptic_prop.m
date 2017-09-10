% proposal distribution from ellipsoid for use with Hermite polynomials
% -----
% [rv,p,p_rat] = herm_elliptic_prop(sample_opt, eval_opt)
% -----
% Input
% -----
% sample_opt = options for sampling inputs
% eval_opt = options for evaluating input and QoI
% ------
% Output
% ------
% rv = sampled inputs
% p = probability of sample in orthogonality distribution. Used for diagnostics
% p_rat = ratio of proposal distribution over orthogonality distribution. Used for MCMC sampler
function [rv,p,p_rat] = herm_elliptic_prop(sample_opt, ~)
    r = sqrt(2).*sqrt(2.*sample_opt.order+1); % Can increased size for slight increase in rejection rate and higher accuracy
    rv = randn(1,sample_opt.r_dim);
    rv = rv/norm(rv,2); % Randomly samples from unit sphere
    x = rand^(1/sample_opt.r_dim).*r; % Determines how far inside along radial axes to move
    rv = x.*rv; % Samples appropriately form interior of ellipsoid
    p = sqrt(2*pi)*exp(-0.5*(rv*rv')); % norm pdf
    p_rat = pi^(0.5*sample_opt.r_dim)/gamma(sample_opt.r_dim/2+1)*prod(r); % Volume of ellipsoid
    p_rat = p_rat/p; % p_rat = ratio of proposal distribution to orthogonality distribution
end
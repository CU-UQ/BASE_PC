% proposal distribution from sphere for use with Hermite polynomials
% -----
% [rv,p,p_rat] = herm_sphere_prop(sample_opt, eval_opt)
% -----
% Input
% -----
% sample_opt = options for sampling inputs
%   sample_opt.order = order (total-order) of PCE used for determining radius
%   sample_opt.r_dim = dimension of input
% eval_opt = options for evaluating input and QoI
% ------
% Output
% ------
% rv = sampled inputs
% p = probability of sample in orthogonality distribution. Used for diagnostics
% p_rat = ratio of proposal distribution over orthogonality distribution. Used for MCMC sampler
function [rv,p,p_rat] = herm_sphere_prop(sample_opt, ~)
    r = sqrt(2)*sqrt(2*max(sample_opt.order)+1); % Can increase size for slight increase in rejection rate and higher accuracy
    rv = randn(1,sample_opt.r_dim);
    rv = rv/norm(rv,2);
    x = rand^(1/sample_opt.r_dim)*r;
    rv = x*rv;
    p = sqrt(2*pi)*exp(-0.5*(rv*rv')); % normal pdf
    p_rat = pi^(0.5*sample_opt.r_dim)/gamma(sample_opt.r_dim/2+1)*r^sample_opt.r_dim/p; % Volume of d-ball over p
end
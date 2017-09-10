% proposal distribution based on orthogonality distribution
% -----
% [rv,p,p_rat] = orth_prop(sample_opt, eval_opt)
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
function [rv,p,p_rat] = orth_prop(~, eval_opt)
    [rv,p] = rv_gen(eval_opt);
    p_rat=1; % ratio is 1 for orthogonality distribution
end

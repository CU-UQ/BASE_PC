% p_rat weight function, useful for draws from proposal distribution that
% are not passed through mcmc sampler
% -----
% w = p_rat_w(lhs,sample_opt) % Actually depends on rv not lhs
% -----
% Input
% -----
% p_rat = ratio of proposal distribution to orthogonality distribution
% ------
% Output
% ------
% w = weight value to be paired with lhs
function w = p_rat_w(~, ~, ~, p_rat, ~)
    w = p_rat;
end

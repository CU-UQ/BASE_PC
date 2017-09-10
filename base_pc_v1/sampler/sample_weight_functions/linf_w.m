% weight function with l_infinity norm, useful for l1-coherence optimal sampling
% -----
% w = linf_w(lhs,sample_opt)
% -----
% Input
% -----
% lhs = vector to associate with weight, associated with a single input
% sample_opt = options for sampling inputs
% ------
% Output
% ------
% w = weight value to be paired with lhs
function w = linf_w(lhs, ~, ~, ~, ~)
    w = 1./norm(lhs,inf);
end

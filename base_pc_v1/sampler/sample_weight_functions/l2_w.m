% weight function with l_2 norm, useful for l2-coherence optimal sampling
% -----
% w = l2_w(lhs,sample_opt)
% -----
% Input
% -----
% lhs = vector to associate with weight, associated with a single input
% sample_opt = options for sampling inputs
% ------
% Output
% ------
% w = weight value to be paired with lhs
function w = l2_w(lhs, ~, ~, ~, ~)
     w = 1./norm(lhs,2);
end

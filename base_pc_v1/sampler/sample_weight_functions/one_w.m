% unit weight function, useful for draws from orthogonality distribution
% -----
% w = one_w(lhs,sample_opt)
% -----
% Input
% -----
% lhs = vector to associate with weight, associated with a single input
% sample_opt = options for sampling inputs
% ------
% Output
% ------
% w = weight value to be paired with lhs
function w = one_w(~, ~, ~, ~, ~)
    w = 1;
end

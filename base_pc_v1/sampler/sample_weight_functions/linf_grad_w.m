% weight function with l_inf norm for use with gradients
% -----
% w = linf_grad_w(lhs,sample_opt)
% -----
% Input
% -----
% lhs = vector to associate with weight, associated with a single input
% sample_opt = options for sampling inputs
% ------
% Output
% ------
% w = weight value to be paired with lhs
function w = linf_grad_w(lhs, ~, ~, ~, sample_opt)
    n_rows = sample_opt.n_grad_dims_plus_one;
    n_samps = size(lhs,1)/n_rows;
    w = zeros(n_samps,1);
    set = 1:n_rows;
    for k = 1:n_samps
        w(k) = 1/sqrt(max(sum(lhs(set,:).^2)));
        set = set+n_rows;
    end
end

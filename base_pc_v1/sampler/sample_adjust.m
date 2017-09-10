% adjusts sample for use with basis
% -----
% sample = sample_adjust(basis, sample, sample_opt, eval_opt)
% -----
% Input
% -----
% basis = basis object
% sample_old = sample object
% sample_opt = options for sampling inputs
% eval_opt = options for evaluating input and QoI
% ------
% Output
% ------
% sample = sample object adjusted to basis
function sample = sample_adjust(basis, sample, sample_opt, eval_opt)
    sample.lhs = basis_eval(basis,sample.rv,eval_opt);
    n_rows = min(eval_opt.grad_dim,basis.n_dim)+1;
    for k = 1:sample.n_samps;
        sample.w(k) = sample_opt.w_handle(sample.lhs((1+(k-1)*n_rows:k*n_rows),:), sample.rv(k,:), sample.p(k), sample.p_rat(k), sample_opt);
        sample.w(k) = sample.w(k)/sample.nc;
    end
    sample.wlhs = apply_weights(sample.w,sample.lhs); % Apply weights
    sample.wrhs = apply_weights(sample.w,sample.rhs); % Apply weights
end
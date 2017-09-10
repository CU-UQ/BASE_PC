% constructs sample object through sampling of orthogonality distribution
% -----
% sample = orth_sampler(basis, sample_opt, eval_opt)
% -----
% Input
% -----
% basis = basis object
% sample_opt = options for sampling inputs
%   sample_opt.n_workers = number of workers in pool
%   sample_opt.initial_size = number of samples to generate
% eval_opt = options for evaluating input and QoI
%   eval_opt.max_dim = maximum dimension of the problem
%   eval_opt.qoi_handle = for rhs evaluations
% ------
% Output
% ------
% sample = sample object

function sample = orth_sampler(sample_opt, eval_opt)
    n_samps = ceil(sample_opt.initial_size); % number of samples initially is predefined as often initial basis warrants different sampling rate
    t_samps = ceil(n_samps/sample_opt.n_workers);
    spmd (sample_opt.n_workers)
        rv = zeros(t_samps,eval_opt.max_dim);
        p = zeros(t_samps,1);
        for k = 1:t_samps
            [rv(k,:),p(k)] = rv_gen(eval_opt);
        end
    end
    sample.rv = cell2mat(rv(:));
    sample.rv = sample.rv(1:n_samps,:);
    sample.p = cell2mat(p(:));
    sample.p = sample.p(1:n_samps,:);
    sample.p_rat = ones(n_samps,1); % Sampling from ortho dist always gives p_rat = 1
    sample.rhs = eval_opt.qoi_handle(sample.rv,eval_opt); % rhs evaluations
    sample.n_samps = n_samps;
end

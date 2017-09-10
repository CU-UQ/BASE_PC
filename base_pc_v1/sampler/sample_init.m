% constructs sample opbject
% -----
% sample = sample_init(basis, sample_opt, eval_opt)
% -----
% Input
% -----
% basis = basis object
% sample_opt = options for sampling inputs
% eval_opt = options for evaluating input and QoI
% ------
% Output
% ------
% sample = sample object

function sample = sample_init(basis, sample_opt, eval_opt)
    n_samps = ceil(sample_opt.initial_size); % Number of samples initially is predefined as often initial basis warrants different sampling rate.
    nc_old = 0;
    sample.rej_con = 0; % Initialize value for rej_con (0 implies chain initialization)
    sample.rej_rate = 0; % Initialize value for rej_rate
    sample.n_tries = 0; % Initialize number of tries. Is used to compute rejection rate.
    w = [];
    while true
        sample  = mcmc_sampler(sample, basis, n_samps, sample_opt, eval_opt);
        w = vertcat(w,sample.w); %#ok<AGROW>
        nc_new = norm(w,2)/size(w,1); % Normalizing constant estimate for working with sampling distribution
        if (nc_new-nc_old)/nc_old < 0.001 % Guages convergence of chain and normalizing constant estimation. Convergence criteria loosens as sampling increases.
            sample.nc = nc_new; % This constant is not often needed to be very accurate to insure quality sampling, if expensive to evaluate.
            break
        end
        nc_old = nc_new;
    end
    sample.w = sample.w/sample.nc; % Adjust weights to .nc
    sample.wlhs = sample.wlhs/sample.nc; % Adjust weights to .nc
    sample.rhs = eval_opt.qoi_handle(sample.rv,eval_opt);%qoi_eval(sample_trans.rv,eval_opt)); % Old rhs plus new evaluations
    sample.wrhs = sample.w.*sample.rhs; % Weights applied
end

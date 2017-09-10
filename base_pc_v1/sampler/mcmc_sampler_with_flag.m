% generates samples via mcmc with potential to flow flags in weight evaluation
% -----
% [sample,flag]  = mcmc_sampler_with_flag(basis, n_samps, samp_opt, eval_opt)
% -----
% Input
% -----
% basis = basis object
% n_samps = number of samples to generate
% sample_opt = options for sampling inputs
% eval_opt = options for evaluating input and QoI
% ------
% Output
% ------
% sample = sample object
% flag = flag returned by weight function to indicate insufficient samples for adaptation
function [sample,flag]  = mcmc_sampler_with_flag(basis, n_samps, samp_opt, eval_opt)
    t_samps = ceil(n_samps/samp_opt.n_workers);
    spmd (samp_opt.n_workers)
        [lhs, rv, w, p, flag] = serial_mcmc_init(basis,t_samps, samp_opt, eval_opt);
    end
    if sum(flag ~= 0) % Return if any flags non-zero
        sample = [];
        return
    end
    n_rows = min(eval_opt.grad_dim,basis.n_dim)+1;
    sample.lhs = cell2mat(lhs(:));
    sample.lhs = sample.lhs(1:(n_rows*n_samps),:);
    sample.rv = cell2mat(rv(:));
    sample.rv = sample.rv(1:n_samps,:);
    sample.p = cell2mat(p(:));
    sample.p = sample.p(1:n_samps);
    sample.w = cell2mat(w(:));
    sample.w = sample.w(1:n_samps);
    sample.n_samps = n_samps;
    set1 = 1:n_rows;
    sample.wlhs = zeros(n_samps*n_rows,basis.n_elems);
    for k = 1:n_samps
        sample.wlhs(set1,:) = sample.w(k)*sample.lhs(set1,:);
        set1 = set1+n_rows;
    end
end

function [lhs, rv, w, p, rej_rate, flag]  = serial_mcmc_init(basis, n_samps, samp_opt, eval_opt)
    flag = [];
    rej_rate = 0; % rejection rate turns serial_red
    n_rows = min(eval_opt.grad_dim,basis.n_dim)+1;
    lhs = zeros(n_samps*n_rows,basis.n_elems);
    w = zeros(n_samps,1);
    p = zeros(n_samps,1);
    rv = zeros(n_samps,eval_opt.max_dim);
    cur_samp.rej_con = 0;
    b = [];
    for k = 1:samp_opt.burn_in % Initial opt.burn_in tunes sampler.
        [cur_samp, rej, flag] = mcmc_try(basis, cur_samp, samp_opt);
        if flag
            return
        end
        rej_rate = rej_rate + rej;
    end
    rej_rate = rej_rate/samp_opt.burn_in;
    n_tries = samp_opt.burn_in;
    serial_red = max(1,ceil(samp_opt.log_col_rate/log(rej_rate)));
    serial_red = min(serial_red,samp_opt.burn_in);
    for k = 1:n_samps % Generate samples to keep
        while true % Remove Collision based duplicates
            rej_s = 0;
            for kkk = 1:serial_red; % Reduces serial correlation
                [b, rej_con, rej, flag] = mcmc_try(basis, b, rej_con, samp_opt);
                if flag
                    return
                end
                rej_s = rej_s + rej;
            end
            rej_rate = rej_rate*n_tries/(n_tries+serial_red)+ rej_s/(n_tries+serial_red); % Adjustment to rejection rate
            n_tries = n_tries+serial_red;
            serial_red_old = serial_red;
            serial_red = max(1,ceil(samp_opt.log_col_rate/log(rej_rate)));
            if(rej_s < serial_red_old)
                break
            end
        end
        lhs((1+(k-1)*n_rows:k*n_rows),:) = b.lhs;
        rv(k,:) = b.rv;
        w(k) = b.w;
        p(k) = b.p;
    end
end

function [cur_samp, rej, flag] = mcmc_try(basis, cur_samp, samp_opt, eval_opt)
    [t_rv,t_p, p_rat] = samp_opt.prop_handle(samp_opt);
    t_lhs = basis_eval(basis,t_rv,eval_opt);    
    old_lhs = basis_eval(sample_opt.old_basis,t_rv,eval_opt);
    [t_w,flag] = sample_opt.w_handle(t_lhs, old_lhs, t_rv, t_p, p_rat, sample_opt); % Many things potentially useful for weight evalution
    tc = p_rat^2*(t_w)^(-2);
    rej = 1; % return 1 if rejected
    if rand < tc/cur_samp.tcp % rejection decision
        cur_samp.lhs = t_lhs;
        cur_samp.w = t_w;
        cur_samp.p = t_p;
        cur_samp.rv = t_rv;
        cur_samp.rej_con = tc;
        rej = 0; % return 0 if accepeted
    end
end

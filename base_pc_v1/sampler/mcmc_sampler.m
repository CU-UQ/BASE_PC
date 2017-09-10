% generates samples via mcmc for sample initialization
% -----
% sample  = mcmc_sampler(sample, basis, n_samps, sample_opt, eval_opt)
% -----
% Input
% -----
% sample = sample object
% basis = basis object
% n_samps = number of samples to generate
% samp_opt = options for sampling inputs
% eval_opt = options for evaluating input and QoI
% ------
% Output
% ------
% sample = sample object
function sample  = mcmc_sampler(sample, basis, n_samps, sample_opt, eval_opt)
    t_samps = ceil(n_samps/sample_opt.n_workers);
    spmd (sample_opt.n_workers)
        [lhs, rv, w, p, p_rat, rej_rate, rej_con, rej_con_history, n_tries] = serial_mcmc_init(basis, t_samps, sample.rej_rate, sample.rej_con, sample.n_tries, sample_opt, eval_opt);
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
    sample.p_rat = cell2mat(p_rat(:));
    sample.p_rat = sample.p_rat(1:n_samps,:);
    sample.rej_con_history = cell2mat(rej_con_history(:));
    sample.rej_con_history = sample.rej_con_history(1:n_samps);

    sample.n_samps = n_samps;
    sample.wlhs = apply_weights(sample.w,sample.lhs);
    sample.rej_rate = sum(cell2mat(rej_rate(:)))/sample_opt.n_workers; % rejection rate estimates
    sample.rej_con = cell2mat(rej_con(randi(sample_opt.n_workers))); % rej_con from random chain (accelerates any future sampling)
    sample.n_tries = sum(cell2mat(n_tries(:))); % Number of tries used to compute rejection rate
end

function [lhs, rv, w, p, p_rat, rej_rate, rej_con, rej_con_history, n_tries]  = serial_mcmc_init(basis, n_samps, rej_rate, rej_con, n_tries, samp_opt, eval_opt)
    n_rows = min(eval_opt.grad_dim,basis.n_dim)+1;
    lhs = zeros(n_samps*n_rows,basis.n_elems);    
    rv = zeros(n_samps,eval_opt.max_dim);
    w = zeros(n_samps,1);
    p = zeros(n_samps,1);
    p_rat = zeros(n_samps,1);
    rej_con_history = zeros(n_samps,1);
    cur_samp = [];
    cur_samp.rej_con = rej_con;
    rej_s = 0;
    for k = 1:samp_opt.burn_in % Initial opt.burn_in tunes sampler. These samples are not kept.
        [cur_samp, rej] = mcmc_try(basis, cur_samp, samp_opt, eval_opt);
        rej_s = rej_s + rej;
    end
    rej_rate = rej_rate*n_tries/(n_tries+samp_opt.burn_in)+ rej_s/(n_tries+samp_opt.burn_in);
    n_tries = n_tries + samp_opt.burn_in;
    serial_red = max(1,ceil(samp_opt.log_col_rate/log(rej_rate)));
    serial_red = min(serial_red,samp_opt.burn_in);
    for k = 1:n_samps % Generate samples to keep
        while true % Remove Collision based duplicates
            rej_s = 0;
            for kkk = 1:serial_red; % Reduces serial correlation
                [cur_samp, rej] = mcmc_try(basis, cur_samp, samp_opt, eval_opt);
                rej_s = rej_s + rej;
            end
            rej_rate = rej_rate*n_tries/(n_tries+serial_red)+ rej_s/(n_tries+serial_red); % Adjustment to rejection rate
            n_tries = n_tries+serial_red;
            serial_red_old = serial_red;
            serial_red = max(1,ceil(samp_opt.log_col_rate/log(rej_rate)));
            if(rej_s < serial_red_old) % We don't keep duplicated samples
                break
            end
        end
        lhs((1+(k-1)*n_rows:k*n_rows),:) = cur_samp.lhs;
        rv(k,:) = cur_samp.rv;
        w(k) = cur_samp.w;
        p(k) = cur_samp.p;
        p_rat(k) = cur_samp.p_rat;
        rej_con_history(k) = cur_samp.rej_con;
        rej_con = cur_samp.rej_con;
    end
end

function [cur_samp, rej] = mcmc_try(basis, cur_samp, sample_opt, eval_opt)
    [t_rv, t_p, p_rat] = sample_opt.prop_handle(sample_opt, eval_opt);
    t_lhs = basis_eval(basis,t_rv,eval_opt);
    t_w = sample_opt.w_handle(t_lhs, t_rv, t_p, p_rat, sample_opt); % Many things potentially useful for weight evalution
    tc = (t_w)^(-2)/p_rat;
    rej = 1; % return 1 if rejected
    if rand < tc/cur_samp.rej_con % rejection decision
        cur_samp.lhs = t_lhs;
        cur_samp.w = t_w;
        cur_samp.p = t_p;
        cur_samp.p_rat = p_rat;
        cur_samp.rv = t_rv;
        cur_samp.rej_con = tc;
        rej = 0; % return 0 if accepeted
    end
end

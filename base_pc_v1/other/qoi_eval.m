% call to evaluate qoi
% -----
% rhs = qoi_eval(input,eval_opt)
% -----
% Input
% -----
% input = point where qoi is evaluated
% eval_opt = options for evaluating input and QoI
% ------
% Output
% ------
% qoi = vector of qoi and derivative information if requested, indexed appropriately
function qoi = qoi_eval(input,eval_opt)
    n_evals = size(input,1);
    if eval_opt.grad
        n_rows = eval_opt.grad_dim+1;
        if eval_opt.parallel_qoi_eval
            qoi = zeros(n_evals,n_rows);
            parfor k = 1:n_evals % Probably easiest to reshape qoi for parallelization
                 qoi(k,:) =  eval_opt.qoi_handle(input(k,:),eval_opt)'; %#ok<PFBNS>
            end
            qoi = qoi(:);
        else
            qoi = zeros(n_rows*n_evals,1);
            set = 1:n_rows;
                for k = 1:n_evals
                    qoi(set) = eval_opt.qoi_handle(input(k,:),eval_opt);
                    set = set + n_rows;
                end
        end
    else
        qoi = zeros(n_evals,1);
        if eval_opt.parallel_qoi_eval
            parfor k = 1:n_evals
                qoi(k) = eval_opt.qoi_handle(input(k,:),eval_opt); %#ok<PFBNS>
            end
        else
            for k = 1:n_evals
                qoi(k) = eval_opt.qoi_handle(input(k,:),eval_opt);
            end
        end
    end
end

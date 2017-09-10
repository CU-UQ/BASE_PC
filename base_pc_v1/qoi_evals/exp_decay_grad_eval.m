% Generates QoI having exponential decay with gradient information
% -----
% output = exp_decay_grad_eval(input, eval_opt)
% -----
% Input
% -----
% input = points where QoI is evaluated. May be multiple rows
% eval_opt = options for evaluating input and QoI
% ------
% Output
% ------
% output = evaluated QoI and gradient information
function output = exp_decay_grad_eval(input,eval_opt)
    output = zeros(eval_opt.grad_dim+1,1);
    output(1) = eval_opt.exp_const;
    for k = 1:eval_opt.max_dim
        output(1) = output(1)-eval_opt.decay_init*input(k)/k^exp_decay_exp;
        if k <= eval_opt.grad_dim
            output(k+1) = -eval_opt.decay_init/k^eval_opt.exp_decay_exp;
        end
    end
    output(1) = exp(output(1));
    output(2:(eval_opt.grad_dim+1)) = output(1)*output(2:(eval_opt.grad_dim+1));
end

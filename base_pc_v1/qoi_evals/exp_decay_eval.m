% Generates QoI having exponential decay
% -----
% output = exp_decay_eval(input, eval_opt)
% -----
% Input
% -----
% input = points where QoI is evaluated. May be multiple rows
% eval_opt = options for evaluating input and QoI
% ------
% Output
% ------
% output = evaluated QoI
function output = exp_decay_eval(input,eval_opt)
    output = eval_opt.exp_const;
    for k = 1:eval_opt.max_dim
        output = output-eval_opt.exp_decay_init.*input(:,k)./k.^eval_opt.exp_decay_exp;
    end
    output = exp(output);
end

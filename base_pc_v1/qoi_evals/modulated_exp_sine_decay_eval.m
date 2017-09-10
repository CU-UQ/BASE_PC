% Generates QoI having exponential decay with non-monotonic decay in dimension
% -----
% output = exp_sine_decay_eval(input, eval_opt)
% -----
% Input
% -----
% input = points where QoI is evaluated. May be multiple rows
% eval_opt = options for evaluating input and QoI
% ------
% Output
% ------
% output = evaluated QoI

function output = modulated_exp_sine_decay_eval(input,eval_opt)
    output = eval_opt.exp_const;
    modulation = 1 - 1./(1+exp(-2*eval_opt.modulation_constant*((input(:,1)-0.5).^2+2*(input(:,2)-0.5).^2-0.25))); % Is small outside circle, is large inside, depending on eval_opt.modulation_constant
    for k = 3:eval_opt.max_dim
        output = output-eval_opt.exp_decay_init.*input(:,k).*sin(k)./k.^eval_opt.exp_decay_exp;
    end
    output = modulation.*exp(output);
end

% Generates QoI evaluating a high order  of sum of inputs
% -----
% output = high_order_eval(input, eval_opt)
% -----
% Input
% -----
% input = points where QoI is evaluated. May be multiple rows
% eval_opt = options for evaluating input and QoI
% ------
% Output
% ------
% output = evaluated QoI

function output = high_order_eval(input,eval_opt)
    output = (input(:,1) + input(:,2) - 1).^eval_opt.power; % Is small outside circle, is large inside, depending on eval_opt.modulation_constant
end

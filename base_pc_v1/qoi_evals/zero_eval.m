% Generates QoI that is always zero
% -----
% output = zero_eval(input, eval_opt)
% -----
% Input
% -----
% input = points where QoI is evaluated. May be multiple rows
% eval_opt = options for evaluating input and QoI
% ------
% Output
% ------
% output = evaluated QoI

function output = zero_eval(input,eval_opt) %#ok<INUSD>
    output  = zeros(size(input,1),1);
end

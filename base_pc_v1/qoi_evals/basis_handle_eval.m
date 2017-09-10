% Generates QoI from predetermined basis and coefficients
% -----
% output = basis_handle_eval(input, eval_opt)
% -----
% Input
% -----
% input = points where QoI is evaluated. May be multiple rows
% eval_opt = options for evaluating input and QoI
% ------
% Output
% ------
% output = evaluated QoI
function output = basis_handle_eval(input, eval_opt)
    output = basis_eval(eval_opt.eval_basis,input,eval_opt)*eval_opt.eval_coefs;
end

% Generates QoI from likelihood_eval (using PC approximation)
% -----
% output = multi_basis_handle_eval(input, eval_opt)
% -----
% Input
% -----
% input = points where QoI is evaluated. May be multiple rows
% eval_opt = options for evaluating input and QoI
% ------
% Output
% ------
% output = evaluated QoI

function output = multi_basis_handle_eval(input,eval_opt)
    output = zeros(1,eval_opt.data_dim); % Row vector
    for k = 1:eval_opt.data_dim
        output(k) =  basis_eval(eval_opt.eval_basis{k},input,eval_opt)*eval_opt.eval_coefs{k}; % assume normal pdf with data mean and variances
    end
end
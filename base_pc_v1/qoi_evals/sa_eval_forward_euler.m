% Generates QoI from surface adsoprtion problem
% -----
% output = sa_eval(input, eval_opt)
% -----
% Input
% -----
% input = points where QoI is evaluated. May be multiple rows
% eval_opt = options for evaluating input and QoI
% ------
% Output
% ------
% output = evaluated QoI

function output = sa_eval_forward_euler(input,eval_opt)
    n_inputs = size(input,1);
    alpha = .1 + exp(10.*input(:,1)); %input(1) normal
    gamma = .001 + .001.*exp(10*input(:,2)); %input(2) normal
    rho = 0.9*ones(n_inputs,1);
    beta = 10;

    % Simple Forward Euler Solver for this
    h = 10^(-5); % Step Size
    tint = 0:h:4;
    n_ints = size(tint,2);
    for k = 1:n_ints
        step = h*(alpha.*(1-rho) - gamma.*rho - beta.*(1-rho).^2.*rho);
        step = max(step,-rho);
        step = min(step,1-rho);
        rho = rho + step;
    end
    output = rho;
end
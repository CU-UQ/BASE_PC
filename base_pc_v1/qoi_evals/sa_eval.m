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

function output = sa_eval(input,eval_opt)
    n_inputs = size(input,1);
    alpha = .1 + exp(10.*input(:,1)); %input(1) normal
    gamma = .001 + .001.*exp(10*input(:,2)); %input(2) normal
    r_0 = 0.9*ones(n_inputs,1);

    % Integration time
    tint = 0:0.01:4;

    % ODE45 options
    opt = odeset('RelTol',1e-5);
    [~,rho] = ode45(@(t,r,a,g) surface_abs(r,alpha,gamma),tint,r_0,opt);
    output = rho(end,:)';
end

function drho = surface_abs(rho,alpha,gamma)
    % Defines the surface absorption forcing
    beta = 10;
    drho = alpha.*(1-rho) - gamma.*rho - beta.*(1-rho).^2.*rho;
end
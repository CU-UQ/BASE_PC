% Generates QoI from Franke function
% -----
% output = franke_eval(input, eval_opt)
% -----
% Input
% -----
% input = points where QoI is evaluated. May be multiple rows
% eval_opt = options for evaluating input and QoI
% ------
% Output
% ------
% output = evaluated QoI

function output = franke_eval(input,eval_opt) %#ok<INUSD>
output = 0.75*exp(-0.25*((9*input(:,1)-2).^2 + (9*input(:,2)-2).^2)) + ...
    0.75*exp(-(9*input(:,1)+1).^2./49 - (9*input(:,2)+1)./10) + ...
    0.5*exp(-0.25*((9*input(:,1)-7).^2 + (9*input(:,2)-3).^2)) - ...
    0.2*exp(-(9*input(:,1)-4).^2 - (9*input(:,2)-7).^2);
end

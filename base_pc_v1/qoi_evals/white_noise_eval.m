% Generates QoI that is white noise
% -----
% output = white_noise_eval(input, eval_opt)
% -----
% Input
% -----
% input = points where QoI is evaluated. May be multiple rows
% eval_opt = options for evaluating input and QoI
% ------
% Output
% ------
% output = evaluated QoI
function output = white_noise_eval(input,eval_opt)
    output  = rand(size(input,1),1)*eval_opt.noise_level;
end

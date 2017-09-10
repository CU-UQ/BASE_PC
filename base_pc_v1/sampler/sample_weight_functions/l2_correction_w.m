% weight function for correction of l2 coherence-optimal sampling
% -----
% [w, t_c, flag] = l2_correction_w(lhs_new, lhs_old, sample_opt)
% -----
% Input
% -----
% lhs_new = vector to associate with weight, associated with a single input
% lhs_old = same input from old basis to be corrected from
% sample_opt = options for sampling inputs
% ------
% Output
% ------
% w = weight value to be paired with lhs_new
% t_c = used for accurate correction sampling
% flag = returned if t_c is zero to indicate sample deficiency
function [w, t_c, flag] = l2_correction_w(lhs_new,lhs_old,~,~,~,sample_opt)
    flag = 0;
    w = 1/norm(lhs_new,2); % value of new weight
    wold = norm(lhs_old,2); % inverse of old weight, but we don't need old weight
    t_c = 1/w^2 - sample_opt.ac*wold^2;
    if t_c < 0 % Correction Sampling is impossible with these parameters
        flag = w^2*wold^2; % 1/o.ac that would have given tc = 0.
    end
end

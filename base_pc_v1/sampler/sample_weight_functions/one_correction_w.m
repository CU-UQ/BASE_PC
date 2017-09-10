% weight function for correction of l2 coherence-optimal sampling
% -----
% [w, t_c, flag] = one_correction_w(lhs_new, lhs_old, sample_opt)
% -----
% Input
% -----
% No input for this function, but must meet params requirement
% ------
% Output
% ------
% w = weight value to be paired with lhs_new
% t_c = used for accurate correction sampling
% flag = returned if t_c is zero to indicate sample deficiency
function [w, t_c, flag] = one_correction_w(~,~,~,~,~,~)
    flag = 0;
    w = 1;
    t_c = 1;
end

% weight function for correction of l1 coherence-optimal sampling
% -----
% [w, t_c, flag] = linf_correction_w(lhs_new, lhs_old, sample_opt)
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

function [w, tc, flag] = linf_correction_w(lhs_new,lhs_old,~,~, ~,sample_opt)
    flag = 0;
    w = 1/norm(lhs_new,inf);
    wold = norm(lhs_old,inf);
    tc = 1/w^2 - sample_opt.ac*wold^2;
    if tc < 0
        flag = w^2*wold^2; % 1/o.ac that would have given tc = 0.
    end
end

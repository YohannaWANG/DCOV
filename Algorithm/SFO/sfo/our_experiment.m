function [A] = our_experiment(Sigma)
%OUR_EXPERIMENT Summary of this function goes here
%   Detailed explanation goes here
V = 1:size(Sigma, 1);
F = sfo_fn_logdet(Sigma, V);

A = sfo_min_norm_point(F, V)
end


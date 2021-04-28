% Computes the log determinant of a covariance matrix
% Author: Arnab Bhattacharyya (arnabb@nus.edu.sg)
%
% function H = sfo_fn_logdet(sigma,set) 
% sigma: Covariance Matrix
% set: the subset of rows
%
% Example: F = sfo_fn_logdet(0.5*eye(3)+0.5*ones(3),1:3);

function F = sfo_fn_logdet(sigma,V, i)
F.sigma = sigma;
F.V = V;
F.exception = i;

F.indsA = [];
F.cholA = [];

F = class(F,'sfo_fn_logdet',sfo_fn);
F = set(F,'current_set',-1);

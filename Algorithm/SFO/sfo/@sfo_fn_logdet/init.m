% Implementation by Arnab Bhattacharyya (arnabb@nus.edu.sg), based on tutorial by Andreas Krause (krausea@gmail.com)
% 
% Example: See sfo_fn.m and the tutorial script for more information

function [F,H] = init(F,sset)
sset = sfo_unique_fast(sset);
if ~isequal(sset,get(F,'current_set'))
    F.cholA = chol(F.sigma([sset, F.exception], [sset, F.exception])+(1e-10)*eye(length(sset)+1));
    F.indsA = sset;
    
    if isempty(sset)
        H = log2(F.sigma(F.exception, F.exception));
    else
        H = sum(log2(diag(F.cholA)));
    end
    F = set(F,'current_val',H,'current_set',sset);
else
    H = get(F,'current_val');
end

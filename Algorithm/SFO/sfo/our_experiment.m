function Abest = our_experiment(Sigma)
%OUR_EXPERIMENT Summary of this function goes here
%   Detailed explanation goes here
V = 1:size(Sigma, 1);
Abest = [];
bestval=0;

for i = V
    Vi = V(V~=i);
    F = sfo_fn_logdet(Sigma, Vi, i);
    A = sfo_min_norm_point(F, Vi);
    if (i==1) || (F(A) < bestval)
        Abest = [A, i];
        bestval = F(A);
    end
end


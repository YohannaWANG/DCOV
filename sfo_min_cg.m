function Abest = sfo_min_cg(Sigma)
dir /home/yohanna/Documents/DCOV/dcov/sfo
%OUR_EXPERIMENT Summary of this function goes here
%   Detailed explanation goes here
V = 1:size(Sigma, 1);
Abest = [];
ibest = 0;
bestval=0;

for i = V
    Vi = V(V~=i);
    F = sfo_fn_logdet(Sigma, Vi, i);
    A = sfo_min_norm_point(F, Vi);
    if (i==1) || (F(A) < bestval)
        Abest = A;
        ibest = i;
        bestval = F(A);
    end
end

Abest = [Abest, ibest];


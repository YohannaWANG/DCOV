% Implementation by Arnab Bhattacharyya (arnabb@nus.edu.sg), based on tutor)ial script by Andreas Krause (krausea@gmail.com)
% 
% Example: See sfo_fn.m and the tutorial script for more information

function newScore = inc(F,A,el)
A = sfo_unique_fast(A);
F = init(F,A);

if sum(A==el)>0
    newScore = get(F,'current_val');
    return;
end

B=[A, F.exception];

sigmaXgA = (F.sigma(el,el)-F.sigma(el,B)*(F.cholA\(F.cholA'\F.sigma(B,el))))/2;

H = sum(log2(sigmaXgA));

newScore = get(F,'current_val')+H;
disp([])

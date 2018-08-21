function mttkrp1 = mttkrp_for_parafac2(YY, K, U, n, PARFOR_FLAG)
% SPARTan entry point computing the MTTKRP specifically for the CP-ALS in
% the PARAFAC2 fitting algorithm.

if (n == 1)
    R = size(U{2},2);
else
    R = size(U{1},2);
end

if (n==1)
    mttkrp1 = mttkrp_mode1(YY, K, U, R, PARFOR_FLAG);
elseif (n==2)
    mttkrp1 =  mttkrp_mode2(YY, K, U, R, PARFOR_FLAG);
elseif (n==3)
    mttkrp1 = mttkrp_mode3(YY, K, U, R, PARFOR_FLAG);
end
    
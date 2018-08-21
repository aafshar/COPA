function mttkrp = mttkrp_mode3(YY, K, U, R, PARFOR_FLAG)
% SPARTan MTTKRP w.r.t. mode-3
u1 = U{1};
u2 = U{2};
mttkrp = zeros(K, R);

if (PARFOR_FLAG)
    parfor k=1: K
        mttkrp(k, :) = dot(u1, YY{k} * u2);
    end
else
    for k=1: K
        mttkrp(k, :) = dot(u1, YY{k} * u2);
    end
end
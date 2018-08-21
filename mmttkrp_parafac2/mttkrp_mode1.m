function mttkrp = mttkrp_mode1(YY, K, U, R, PARFOR_FLAG)
% SPARTan MTTKRP w.r.t. mode-1

u2 = U{2}; 
u3 = U{3};
mttkrp = zeros(R, R);

if (PARFOR_FLAG)
    parfor k=1: K
        mttkrp = mttkrp + bsxfun(@times, u3(k, :), YY{k} * u2);
    end
else
    for k=1: K
        mttkrp = mttkrp + bsxfun(@times, u3(k, :), YY{k} * u2);
    end
end

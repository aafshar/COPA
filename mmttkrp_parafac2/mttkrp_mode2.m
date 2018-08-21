function mttkrp = mttkrp_mode2(YY, K, U, R, PARFOR_FLAG)
% SPARTan MTTKRP w.r.t. mode-2
u1 = U{1}; 
u3 = U{3};
I = size(U{2}, 1);
mttkrp = zeros(I, R);

if (PARFOR_FLAG)
    parfor k=1: K
        un = unique_col_ind(YY{k});
        partial_sum = zeros(I, R);
        partial_sum(un, :) = bsxfun(@times, u3(k, :), YY{k}(:, un)' * u1);    
        mttkrp = mttkrp + partial_sum;
    end
else
    for k=1: K
        un = unique_col_ind(YY{k});
        partial_sum = zeros(I, R);
        partial_sum(un, :) = bsxfun(@times, u3(k, :), YY{k}(:, un)' * u1);    
        mttkrp = mttkrp + partial_sum;
    end
end




%%
% u1 = U{1}; 
% u3 = U{3};
% mttkrp = zeros(size(U{2}, 1), R);
% parfor k=1: K
%     mttkrp = mttkrp + bsxfun(@times, u3(k, :), YY{k}' * u1);
% end
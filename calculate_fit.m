function [ fit ] = calculate_fit( X,Q,H,W,V,normX,K,PARFOR_FLAG )
%Calculate fit for parafac2 problem
    fit=0;
    if (PARFOR_FLAG)
        parfor k = 1:K
         M   = (Q{k}*H)*diag(W(k,:))*V';
         fit = fit + sum(sum( (X{k} - M ).^2));
        end
        fit=1-(fit/normX);
    else
        for k = 1:K
         M   = (Q{k}*H)*diag(W(k,:))*V';
         fit = fit + sum(sum( (X{k} - M ).^2));
        end
        fit=1-(fit/normX);

    end

end


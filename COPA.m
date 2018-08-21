function [fit,FIT_TIME,U,V,W] = COPA( X,R,conv_tol,seed,PARFOR_FLAG,normX,Constraints )
% This code contains the core part of COPA framework 

tStart=tic;
FIT_TIME=[];
ROOTPATH = '';

J=size(X{1}, 2); %  number of features (variables)
K = max(size(X));% number of subjects
Q=cell(K,1);
YY=cell(K,1);

rng(seed); % initilizing the modes based on some seed
V = rand(J,R);
W = rand(K,R);
H = eye(R);
prev_fit=0; fit=1;
itr=0;


while(abs(fit-prev_fit)>conv_tol*prev_fit)
    itr=itr+1;
    if (PARFOR_FLAG)
        parfor k=1:K
            Qk = H*diag(W(k,:))*(X{k}*V)';%(A'*X{k}');
            Q{k} = Qk'*psqrt(Qk*Qk');
            YY{k} = sparse(Q{k}'*X{k}); 
        end
    else
        for k=1:K

            Qk = H*diag(W(k,:))*(X{k}*V)';%(A'*X{k}');
            Q{k} = Qk'*psqrt(Qk*Qk');
            YY{k} = sparse(Q{k}'*X{k});
            
            
        end
    end
    t_ten=tic;
    [ Tensor ] = COPA_optimizer( YY, R,  'maxiters', 1 ,'init', {H, V, W},'Constraints',Constraints,'PARFOR_FLAG',PARFOR_FLAG );
    toc(t_ten)
    H=Tensor{1};
    V=Tensor{2};
    W=Tensor{3};

    prev_fit=fit;
    fit=calculate_fit(X,Q,H,W,V,normX,K,PARFOR_FLAG)
        
       
        

    tEnd = toc(tStart);
    FIT_TIME(itr,1)=tEnd;
    FIT_TIME(itr,2)=fit;
    

    
end
%construct the U_k
U=cell(K,1);
if(PARFOR_FLAG)
    parfor k=1:K
        U{k}=Q{k}*H;
    end
else
    for k=1:K
        U{k}=Q{k}*H;
    end
end


%plot U_k for a random subject
subject_number=1;
for r=1:R
    plot([1:size(U{subject_number},1)],U{subject_number}(:,r))
    hold on;
end
ylabel("Value")
xlabel("Number of observations")
title("plot U{k} where k is 1")
end

function X = psqrt(A,tol)
% Produces A^(-.5) even if rank-problems

[U,S,V] = svd(A,0);
if min(size(S)) == 1
    S = S(1);
else
    S = diag(S);
end
if (nargin == 1)
    tol = max(size(A)) * S(1) * eps;
end
r = sum(S > tol);
if (r == 0)
    X = zeros(size(A'));
else
    S = diag(ones(r,1)./sqrt(S(1:r)));
    X = V(:,1:r)*S*U(:,1:r)';
end

end



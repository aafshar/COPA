function [ fit,FIT_TIME,U,V,W ] = Smooth_COPA( X,R,conv_tol,seed,PARFOR_FLAG,normX,Constraints,GAP )
%Implementation of smooth parafac 2 where smoothness impose on mode U_k
% If GAP=1 then the smoothness considers the gap between two time stamps

tStart=tic;
FIT_TIME=[];

J=size(X{1}, 2); %  number of features
K = max(size(X));% number of subjects
if(GAP==1)
    %fid = fopen('gaps_per_subject.csv','rt');
    % the file format should be like this:
    %each line is related to a subject. (4,7,12,19,39,45,......)
    %each number is a time stamp.
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
F={};%containts all basis functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Ra=cell(K,1);
num_of_basis_fun=7;
if(PARFOR_FLAG)
    parfor k=1:K
        %create B spline for each subject
        if(GAP==1) % Incorporating the observational gap for creating the smooth function.
            knots = str2num(fgetl(fid)) %read the next line of the days of a patient
            knots=knots-(knots(1)-1);
            patient_dist=knots/knots(end);
            knots=knots*size(X{k},1);
            
        else
            knots=1:size(X{k},1);
        end
        
        F{k}=MSplineBasis([knots], num_of_basis_fun,3, [knots(1) knots(end)] ); 
        [u,~,~]=svd(F{k},'econ');
        Ra{k}=u*u';

    end
else
    for k=1:K
        if(GAP==1)
            knots = str2num(fgetl(fid)); %read the next line of the days of a patient
            knots=knots-(knots(1)-1);
            knots=knots/knots(end);
            knots=knots*size(X{k},1);
            
        else
            knots=1:size(X{k},1);
        end
        F{k}=MSplineBasis([knots], num_of_basis_fun,3, [knots(1) knots(end)] ); 
        [u,~,~]=svd(F{k},'econ');
        Ra{k}=u*u';

    end
end

rng(seed)
H=rand(R,R);
V=rand(J,R);
W=rand(K,R);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
itr=0;
Rknew=cell(K,1);
Xtilde=cell(K,1);
prev_fit=0; fit=1;


while(abs(fit-prev_fit)>conv_tol*prev_fit)
    prev_fit=fit;
    if(PARFOR_FLAG)
        parfor k=1:K

            [u,~,v]=svd(Ra{k}*X{k}*V*(diag(W(k,:)))*H','econ');
            Rknew{k}=u*v';
            Xtilde{k} = sparse(Rknew{k}'*X{k});

        end
    else
        for k=1:K

            [u,~,v]=svd(Ra{k}*X{k}*V*(diag(W(k,:)))*H','econ');
            Rknew{k}=u*v';
            Xtilde{k} = sparse(Rknew{k}'*X{k});

        end
    end


        
    [ Tensor ]  = COPA_optimizer( Xtilde, R,  'maxiters', 1 ,'init', {H, V, W},'Constraints',Constraints,'PARFOR_FLAG',PARFOR_FLAG );
    H=Tensor{1};
    V=Tensor{2};
    W=Tensor{3};
     
    itr=itr+1;
    
    fit=calculate_fit(X,Rknew,H,W,V,normX,K,PARFOR_FLAG)
    tEnd = toc(tStart);
    FIT_TIME(itr,1)=tEnd;
    FIT_TIME(itr,2)=fit;

end


U=cell(K,1);
if(PARFOR_FLAG)
    parfor k=1:K
        U{k}=Rknew{k}*H;
    end
else
    for k=1:K
        U{k}=Rknew{k}*H;
    end
end

%plot U_k for a random subject w/o gap.
subject_number=1;
if(GAP==1)
    M = csvread('gaps_per_subject.csv');
    temp=M(subject_number,:);
    temp(temp==0) = [];
    for r=1:R
      plot(temp,U{subject_number}(:,r))
      hold on;
    end

else
    for r=1:R
      plot([1:size(U{subject_number},1)],U{subject_number}(:,r))
      hold on;
    end
end
ylabel("Value")
xlabel("Number of observations")
title("plot U{k} where k is 1")



end

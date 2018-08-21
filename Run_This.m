clc;
clear;

%Import neccessary packages
addpath(genpath('./N-Way'));
addpath(genpath('./tensor_toolbox'));
addpath(genpath('./mmttkrp_parafac2'));


load('parafac2_problem.mat', 'X'); % load the data set
K = max(size(X)); %number of subjects (K)
R = 5; %number of factors or components
PARFOR_FLAG=0;
normX=claculate_norm(X,K,PARFOR_FLAG); %Calculate the norm of the input X
conv_tol=1e-4; %converegance tolerance
Smoothness=1; %if smoothness is 1, the U_k will be smooth otherwise not.

%Defining the constraints here:
%you can use the folowing constraints in COPA:
%1-non-negativity on H, S_k, and V
%2-Smoothness on U_k
%3-sparsity constraint (l1 and l0) on H, S_k, V

%the first element impose a constraint on H
%the second element impose a constraint on V
%the third element impose a constraint on W (S_k)
Constraints={'nonnegative', 'nonnegative','nonnegative'};
seed=2;
if Smoothness==1
    tStart=tic;
    GAP=0;
    [fit,FIT_TIME]=Smooth_COPA(X,R,conv_tol,seed,PARFOR_FLAG,normX,Constraints,GAP );
    tEnd = toc(tStart);
else
    tStart=tic;
    [fit,FIT_TIME]=COPA(X,R,conv_tol,seed,PARFOR_FLAG,normX,Constraints);
    tEnd = toc(tStart);
end



%plot the FIT-Time figure
figure;
plot(FIT_TIME(:,1),FIT_TIME(:,2));
xlabel("TIME")
ylabel("FIT")




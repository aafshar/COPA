function [ H] = COPA_optimizer( X, R, varargin )
%COPA framework for optimizing H,W (S_k) and V.

K=size(X,1);%number of subjects
size_ten = [size(X{1}, 1), size(X{1}, 2), K]; %Tensor Size
N=length(size_ten);

%% Set algorithm parameters from input or by using defaults
params = inputParser;
params.addParamValue('maxiters',50,@(x) isscalar(x) & x > 0);
params.addParamValue('init', 'random', @(x) (iscell(x) || ismember(x,{'random','nvecs'})));
params.addParamValue('Constraints',{'nonnegative','nonnegative','nonnegative' },@(x) length(x)==3);
params.addParamValue('PARFOR_FLAG', 0, @isscalar);
params.parse(varargin{:});

%% Copy from params object
maxiters = params.Results.maxiters;
init = params.Results.init;
Constraints = params.Results.Constraints;
PARFOR_FLAG=params.Results.PARFOR_FLAG;
GG=cell(N,1);
U=cell(N,1);
if (PARFOR_FLAG)
    parfor n = 1:N
        GG{n} = init{n}'*init{n};
        U{n}=zeros(size_ten(n),R);
    end
else
    for n = 1:N
        GG{n} = init{n}'*init{n};
        U{n}=zeros(size_ten(n),R);
    end
end


for itr = 1:1
    for n =1:N
       [init U, GG]=fastADMM( X, init, U, n, GG, Constraints,PARFOR_FLAG);

    end
end
H=init;


end

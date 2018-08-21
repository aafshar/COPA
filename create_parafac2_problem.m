function [X, Mat] = create_parafac2_problem(K, J, F, sparsity, MAXSAMPLES, PARFOR_FLAG)
rng('default');

C = rand(K, F);
A = rand(J, F);
L=rand(F,50);
Mat=C*L;
H=orth(orth(rand(F))');
H(H<0);
assert(rank(H)==size(H, 1));

P=cell(1, K);
if (PARFOR_FLAG)
    parfor i=1: K
        P{i}=orth(rand(MAXSAMPLES, F));
    end
else
    for i=1: K
        P{i}=orth(rand(MAXSAMPLES, F));
    end
end

X = cell(1, K);
totalnnz = 0;
if (PARFOR_FLAG)
    parfor i=1: K
        X{i}=(A*diag(C(i,:))*(P{i}*H)')' ;
        %sparsifier = sprand(MAXSAMPLES, J, sparsity);
        %sparsifier = (sparsifier>0);
        %X{i}  = X{i}.*sparsifier;
        %X{i} = sparse(X{i});
        X{i} = X{i}(any(X{i}, 2), :);
        totalnnz = totalnnz + nnz(X{i});
    end
else
    for i=1: K
        X{i}=(A*diag(C(i,:))*(P{i}*H)')' ;
        %sparsifier = sprand(MAXSAMPLES, J, sparsity);
        %sparsifier = (sparsifier>0);
        %X{i}  = X{i}.*sparsifier;
        %X{i} = sparse(X{i});
        X{i} = X{i}(any(X{i}, 2), :);
        totalnnz = totalnnz + nnz(X{i});
    end
end
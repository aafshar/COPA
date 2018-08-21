function [ H, U, GG, itr ] = fastADMM( Y, H, U, d, GG, Constraints,PARFOR_FLAG)
% ADMM iterates to solve one mode of tensor Y. 
%based on formula 8 in the paper.

K=size(Y,1);
itr=0;
[ ~, k ] = size(H{d});
G = ones(k,k); prod = [ 1:d-1, d+1:length(GG) ];
for dd = prod
    G = G .* GG{dd}; 
end

rho =min( 1e-3,trace(G)/k);
L = chol( G + (rho)*eye(k), 'lower' );

%F = mttkrp( Y, H, d );
F=mttkrp_for_parafac2(Y, K, H, d, PARFOR_FLAG);
tol = 1e-2;
Hd = H{d}; Ud = U{d};
for itr = 1:1 %you can change the number of inner iterations
    H0 = Hd;

    Ht   = L'\ ( L\ ( F + rho*(Hd+Ud) )' );
    Hd = proxr( Ht'-Ud, Constraints, d, rho);
    Ud = Ud + Hd - Ht';

end

U{d} = Ud;

    
H{d} = Hd;   
GG{d} = Hd'*Hd;
end


function H = proxr( Hb, Constraints, d, rho )
    switch Constraints{d}
        case 'nonnegative'
            H = max( 0, Hb );
        case 'l1'
            l1_regul=0.0000085;
            H = sign( Hb ) .* max( 0, abs(Hb) - (l1_regul/rho) ); 
        case 'l0'
            l0_regul=5;
            Hb(Hb <= l0_regul) = 0; 
            H=Hb;

    end
end

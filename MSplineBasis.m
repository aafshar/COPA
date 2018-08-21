function [ Xmat ] = MSplineBasis( x, df,degree,boundry )
%This function produce basis functions for a subject
%x is the input. 
% degree is the degree of freedom for Bspline function
x(end)=x(end)-1e-5;
order=degree+1;
nIknots = df - order;
knots=linspace(0,1,nIknots+2);
knots=knots(2:end-1);
knots=quantile(x,knots);
Aknots=sort([repmat(boundry,1,order) knots]);
n=length(x);
nk=length(Aknots);
knotdiferrence=Aknots(2:nk)-Aknots(1:(nk - 1));
knotdiferrence(knotdiferrence==0)=1;
Xmat=(bsxfun( @ge, x', Aknots(1:(nk-1)) ) & bsxfun( @lt, x',Aknots(2:nk) )).*(repmat((1./knotdiferrence),n,1));

for k=2:order
    X1new=bsxfun( @minus, x', Aknots(1:(nk-k))).*Xmat(:,1:nk-k);
    X2new=bsxfun( @plus, -x', Aknots((1+k):nk)).*Xmat(:,2:nk-k+1);
    knotdiferrence= Aknots((1 + k):nk) - Aknots(1:(nk - k));
    knotdiferrence(knotdiferrence==0)=1;
    Xmat=(k/(k - 1)) * (X1new + X2new) .*repmat((1./knotdiferrence(1:nk-k)),n,1);
end


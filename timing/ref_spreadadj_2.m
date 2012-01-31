function cadj=ref_spreadadj_2(coef)
%REF_SPREADADJ_2  Symbol of adjoint spreading function.
%   Usage: cadj=ref_spreadadj_2(c,number);
%
%   Development version by FJ for comparison of different implementations
%   cadj=SPREADADJ(c) will compute the symbol cadj of the spreading 
%   operator that is the adjoint of the spreading operator with symbol c. 
%
%   This is a simple implementation of the direct formula
%

L=size(coef,1);

% It is possible to find a direct one to one mapping between the 
% coefficients coef and cadj which can be written in the following 
% simple and readable form :
        
cadj=zeros(L);
for ii=0:L-1
  for jj=0:L-1
    cadj(ii+1,jj+1)=conj(coef(mod(L-ii,L)+1,mod(L-jj,L)+1))...
        *exp(-i*2*pi*ii*jj/L);
  end;
end;
        



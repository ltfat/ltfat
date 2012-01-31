function cadj=ref_spreadadj_1(coef);
%REF_SPREADADJ  Symbol of adjoint preading function.
%   Usage: cadj=ref_spreadadj_1(c);
%
%   cadj=SPREADADJ(c) will compute the symbol cadj of the spreading
%   operator that is the adjoint of the spreading operator with symbol c.
%
%   The algorithm uses an explit algorithm to compute the entries.
  

L=size(coef,1);

cadj=zeros(L);
for ii=0:L-1
  for jj=0:L-1
    cadj(ii+1,jj+1)=conj(coef(mod(L-ii,L)+1,mod(L-jj,L)+1))...
        *exp(-i*2*pi*ii*jj/L);
  end;
end;





function cadj=ref_spreadadj_6(coef)
%REF_SPREADADJ_6  Symbol of adjoint spreading function.
%   Usage: cadj=ref_spreadadj_6(c,number);
%
%   Development version by FJ for comparison of different implementations
%   cadj=SPREADADJ(c) will compute the symbol cadj of the spreading 
%   operator that is the adjoint of the spreading operator with symbol c. 
%
%   Implementation for sparse matrix without loop (faster)

L=size(coef,1);
  
[row,col,val]=find(coef);
        
% Optimization note : see the corresponding note in case 5
temp=exp((-i*2*pi/L)*(0:L-1)');

ii=mod(L-row+1, L);
jj=mod(L-col+1, L);
cadj=sparse(ii+1,jj+1,conj(val).*temp(mod(ii.*jj,L)+1),L,L);        




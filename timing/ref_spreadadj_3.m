function cadj=ref_spreadadj_3(coef)
%REF_SPREADADJ_3  Symbol of adjoint spreading function.
%   Usage: cadj=ref_spreadadj_3(c,number);
%
%   Development version by FJ for comparison of different implementations
%   cadj=SPREADADJ(c) will compute the symbol cadj of the spreading 
%   operator that is the adjoint of the spreading operator with symbol c. 
%
%   This is an improved implementation of the direct formula with low memory
%   needs
%
%   This implementation uses the direct formula given in case 2 with 
%   the following Optimizations :
%
%-     Avoiding mod : In the loop of case 2, we see that 
%      mod(L-ii,L)~=L-ii only for ii==0 (idem for jj), so we can
%      remove the mod by processing separetly the cases ii==0 or
%      jj==0.
%
%-     Precomputation of exp : In the loop of case 2, we see that we
%      compute many time complex exponential terms with the same 
%      values. Using precomputation and modulo, we can reduce the
%      computation time
%-
%

L=size(coef,1);

cadj=zeros(L);

% Proceesing for ii==0 or jj==0
cadj(1,1)=conj(coef(1,1));
cadj(2:end,1)=conj(coef(end:-1:2,1));
cadj(1,2:end,1)=conj(coef(1,end:-1:2));

% Proceesing for ii~=0 and jj~=0

% Precomputation for exponential term
% Optimization note : here we are computing the Lth root of unity 
% which have many known special properties and symetries that 
% could be exploited to highly reduce the computation for the 
% following line (this is not a critical part but some 
% improvements Optimizations are easy to identify)
temp=exp((-i*2*pi/L)*(0:L-1));

for ii=1:L-1
  for jj=1:L-1
    cadj(ii+1,jj+1)=conj(coef(L-ii+1,L-jj+1))...
        *temp(mod(ii*jj,L)+1);
  end;
end;



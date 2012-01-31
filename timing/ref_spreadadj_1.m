function cadj=ref_spreadadj_1(coef)
%REF_SPREADADJ_1  Symbol of adjoint spreading function.
%   Usage: cadj=ref_spreadadj_1(c,number);
%
%   Development version by FJ for comparison of different implementations
%   cadj=SPREADADJ(c) will compute the symbol cadj of the spreading 
%   operator that is the adjoint of the spreading operator with symbol c. 
%
%   This is the implementation previously in the toolbox, 
%   with the addition of an initialisation of cadj before the loop 
%   leading to a huge speed improvement :

L=size(coef,1);

% Matlab cannot handle an FFT of a sparse matrix.
coef=ifft(full(coef));

cadj = zeros(L);
for ii=0:L-1
  for jj=0:L-1
    cadj(ii+1,jj+1)=conj(coef(mod(ii-jj,L)+1,mod(-jj,L)+1));
  end;
end;

cadj=fft(full(cadj));



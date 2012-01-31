function coef=ref_spreadfun(T)
%REF_SPREADFUN  Spreading function.
%   Usage:  c=ref_spreadfun(T);
%
%   REF_SPREADFUN(T) computes the spreading function of the operator T. The
%   spreading function represent the operator T as a weighted sum of
%   time-frequency shifts. See the help text for SPREADOP for the exact
%   definition.
%
%   SEE ALSO:  SPREADOP, TCONV, SPREADINV, SPREADADJ

L=size(T,1);

coef=zeros(L);
for ii=0:L-1
  for jj=0:L-1
    coef(ii+1,jj+1)=T(ii+1,mod(ii-jj,L)+1);
  end;
end;

coef=fft(coef)/L;




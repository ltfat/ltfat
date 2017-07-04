function [V,D]=spreadeigs(K,coef);
%SPREADEIGS  Eigenpairs of Spreading operator
%   Usage: h=spreadeigs(K,c);
%
%   `spreadeigs(K,c)` computes the *K* largest eigenvalues and eigen-
%   vectors of the spreading operator with symbol *c*.
%
%   See also:  tconv, spreadfun, spreadinv, spreadadj

complainif_argnonotinrange(nargin,2,2,mfilename);

if ndims(coef)>2 || size(coef,1)~=size(coef,2)
    error('Input symbol coef must be a square matrix.');
end;

L=size(coef,1);

% This version explicitly constucts the matrix representation T
% and then applies this matrix as the final step.
coef=fft(coef);
  
T=zeros(L);
for ii=0:L-1
  for jj=0:L-1
    T(ii+1,jj+1)=coef(ii+1,mod(ii-jj,L)+1);
  end;
end;

if nargout==2
  doV=1;
else
  doV=0;
end;

if doV
  [V,D]=eig(T);
else
  D=eig(T);
end;

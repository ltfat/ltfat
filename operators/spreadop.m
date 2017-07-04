function h=spreadop(f,coef)
%SPREADOP  Spreading operator
%   Usage: h=spreadop(f,c);
%
%   `spreadop(f,c)` applies the operator with spreading function *c* to the
%   input *f*. *c* must be square.
%
%   `spreadop(f,c)` computes the following for *c* of size $L \times L$:
% 
%   ..         L-1 L-1 
%     h(l+1) = sum sum c(m+1,n+1)*exp(2*pi*i*l*m/L)*f(l-n+1)
%              n=0 m=0
%
%   .. math:: h\left(l+1\right)=\sum_{n=0}^{L-1}\sum_{m=0}^{L-1}c\left(m+1,n+1\right)e^{2{\pi}ilm/L}f\left(l-n+1\right)
%
%   where $l=0,\ldots,L-1$ and $l-n$ is computed modulo *L*.
%
%   The combined symbol of two spreading operators can be found by using
%   `tconv`. Consider two symbols *c1* and *c2* and define *f1* and *f2* by::
%
%     h  = tconv(c1,c2)
%     f1 = spreadop(spreadop(f,c2),c1);
%     f2 = spreadop(f,h);
%
%   then *f1* and *f2* are equal.
%
%   See also:  tconv, spreadfun, spreadinv, spreadadj
%
%   References: feko98

%   AUTHOR: Peter L. Søndergaard
%   TESTING: TEST_SPREAD
%   REFERENCE: REF_SPREADOP

complainif_argnonotinrange(nargin,2,2,mfilename);

if ndims(coef)>2 || size(coef,1)~=size(coef,2)
    error('Input symbol coef must be a square matrix.');
end;

L=size(coef,1);

% Change f to correct shape.
[f,Ls,W,wasrow,remembershape]=comp_sigreshape_pre(f,'DGT',0);

f=postpad(f,L);

h=zeros(L,W);

if issparse(coef) && nnz(coef)<L

    % The symbol is so sparse that the straighforward definition is
    % the fastest way to apply it.

    [mr,nr,cv]=find(coef);
    h=zeros(L,W);

    % We need mr and nr to be zero-indexed
    mr=mr-1;
    nr=nr-1;

    % This is the basic idea of the routine below
    %for ii=1:length(mr)
    %  for l=0:L-1
    %	h(l+1,:)=h(l+1,:)+cv(ii)*exp(-2*pi*i*l*mr(ii)/L)*f(mod(l-nr(ii),L)+1,:);
    %  end;
    %end;

    l=(-2*pi*i*(0:L-1)/L).';
    for ii=1:length(mr)
        bigmod=repmat(exp(l*mr(ii)),1,W);
        h=h+cv(ii)*(bigmod.*circshift(f,nr(ii)));
    end;

else

  % This version only touches coef one column at a time, and it suited
  % if coef is sparse.
  
  if issparse(coef)
    for n=0:L-1
      % The 'full' below is required for Matlab compatibility, as
      % Matlab refuses to do an ifft of a sparse matrix.
      cf=ifft(full(coef(:,n+1)))*L;
      modind=mod((0:L-1).'-n,L)+1;
      h=h+repmat(cf,1,W).*f(modind,:);
    end;            
  else
    for n=0:L-1
      cf=ifft(coef(:,n+1))*L;
      modind=mod((0:L-1).'-n,L)+1;
      h=h+repmat(cf,1,W).*f(modind,:);
    end;                            
  end;
  
end;


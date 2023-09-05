function fd=pderiv(f,dim,difforder)
%PDERIV   Derivative of smooth periodic function
%   Usage:  fd=pderiv(f);
%           fd=pderiv(f,dim);
%           fd=pderiv(f,dim,difforder);
%
%   `pderiv(f)` will compute the derivative of *f* using a 4th order
%   centered finite difference scheme. *f* must have been obtained by a
%   regular sampling. If *f* is a matrix, the derivative along the columns
%   will be found.
%
%   `pderiv(f,dim)` will do the same along dimension *dim*.
%
%   `pderiv(f,dim,difforder)` uses a centered finite difference scheme of
%   order *difforder* instead of the default.
%
%   `pderiv(f,dim,Inf)` will compute the spectral derivative using a DFT.
%
%   `pderiv` assumes that *f* is a regular sampling of a function on the
%   torus $[0,1)$. The derivative of a function on a general torus $[0,T)$
%   can be found by scaling the output by $1/T$. 

% Assert correct input.

complainif_argnonotinrange(nargin,1,3,mfilename);

if nargin==1
  dim=[];
end;

if nargin<3
  difforder=4;
end;

[f,L,Ls,W,dim,permutedsize,order]=assert_sigreshape_pre(f,[],dim,'PDERIV');

switch(difforder)
 case 2
  fd = L*(circshift(f,-1)-circshift(f,1))/2;
 case 4
  fd = L*(-circshift(f,-2)+8*circshift(f,-1)-8*circshift(f,1)+ ...
          circshift(f,2))/12;
 case Inf
  n=fftindex(L,0);
  n=repmat(n,1,W);
  
  fd=2*pi*ifft(i*n.*fft(f));

  if isreal(f)
    fd=real(fd);
  end;

 otherwise
  error('The specified differentation order is not implemented.');
end;

fd=assert_sigreshape_post(fd,dim,permutedsize,order);


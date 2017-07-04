function c=dstiv(f,L,dim)
%DSTIV  Discrete Sine Transform type IV
%   Usage:  c=dstiv(f);
%           c=dstiv(f,L);
%           c=dstiv(f,[],dim);
%           c=dstiv(f,L,dim);
%
%   `dstiv(f)` computes the discrete sine transform of type IV of the input
%   signal *f*. If *f* is a matrix, then the transformation is applied to
%   each column. For N-D arrays, the transformation is applied to the first
%   non-singleton dimension.
%
%   `dstiv(f,L)` zero-pads or truncates f to length *L* before doing the
%   transformation.
%
%   `dstiv(f,[],dim)` applies the transformation along dimension *dim*. 
%   `dstiv(f,L,dim)` does the same, but pads or truncates to length *L*.
%   
%   The transform is real (output is real if input is real) and
%   it is orthonormal. It is its own inverse.
%
%   Let *f* be a signal of length *L* and let `c=dstiv(f)`. Then
%
%   ..                     L-1
%     c(n+1) = sqrt(2/L) * sum f(m+1)*sin(pi*(n+.5)*(m+.5)/L) 
%                          m=0 
%
%   .. math:: c\left(n+1\right)=\sqrt{\frac{2}{L}}\sum_{m=0}^{L-1}f\left(m+1\right)\sin\left(\frac{\pi}{L}\left(n+\frac{1}{2}\right)\left(m+\frac{1}{2}\right)\right)
%
%   Examples:
%   ---------
%
%   The following figures show the first 4 basis functions of the DSTIV of
%   length 20:::
%
%     % The dstiv is its own adjoint.
%     F=dstiv(eye(20));
%
%     for ii=1:4
%       subplot(4,1,ii);
%       stem(F(:,ii));
%     end;
%
%   See also:  dstii, dstiii, dctii
%
%   References: rayi90 wi94

%   AUTHOR: Peter L. Søndergaard
%   TESTING: TEST_PUREFREQ
%   REFERENCE: REF_DSTIV

complainif_argnonotinrange(nargin,1,3,mfilename);

if nargin<3
  dim=[];
end;

if nargin<2
  L=[];
end;
    
[f,L,Ls,W,dim,permutedsize,order]=assert_sigreshape_pre(f,L,dim,'DSTIV');

if ~isempty(L)
  f=postpad(f,L);
end;

c = comp_dst(f,4);

c=assert_sigreshape_post(c,dim,permutedsize,order);

% This is a slow, but convenient way of expressing the algorithm.
%R=1/sqrt(2)*[diag(exp(-(0:L-1)*pi*i/(2*L)));...
%	     flipud(diag(-exp((1:L)*pi*i/(2*L))))];

%c=i*(exp(-pi*i/(4*L))*R.'*fft(R*f)/sqrt(2*L));


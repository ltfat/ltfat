function c=dstii(f,L,dim)
%DSTII  Discrete Sine Transform type II
%   Usage:  c=dstii(f);
%           c=dstii(f,L);
%           c=dstii(f,[],dim);
%           c=dstii(f,L,dim);
%
%   `dstii(f)` computes the discrete sine transform of type II of the
%   input signal *f*. If *f* is multi-dimensional, the transformation is
%   applied along the first non-singleton dimension.
%
%   `dstii(f,L)` zero-pads or truncates *f* to length *L* before doing the
%   transformation.
%
%   `dstii(f,[],dim)` or `dstii(f,L,dim)` applies the transformation along
%   dimension *dim*.
%
%   The transform is real (output is real if input is real) and orthonormal.
%
%   The inverse transform of |dstii| is |dstiii|.
%
%   Let *f* be a signal of length *L*, let `c=dstii(f)` and define the vector
%   *w* of length *L* by
%
%   .. w = [1 1 1 1 ... 1/sqrt(2)]
%
%   .. math:: w\left(n\right)=\begin{cases}\frac{1}{\sqrt{2}} & \text{if }n=L-1\\1 & \text{otherwise}\end{cases}
%
%   Then 
%
%   ..                     L-1
%     c(n+1) = sqrt(2/L) * sum w(n+1)*f(m+1)*sin(pi*n*(m+.5)/L) 
%                          m=0 
%
%   .. math:: c\left(n+1\right)=\sqrt{\frac{2}{L}}\sum_{m=0}^{L-1}w\left(n\right)f\left(m+1\right)\sin\left(\frac{\pi}{L}n\left(m+\frac{1}{2}\right)\right)
%
%   Examples:
%   ---------
%
%   The following figures show the first 4 basis functions of the DSTII of
%   length 20:::
%
%     % The dstiii is the adjoint of dstii.
%     F=dstiii(eye(20));
%
%     for ii=1:4
%       subplot(4,1,ii);
%       stem(F(:,ii));
%     end;
%
%   See also:  dctii, dstiii, dstiv
%
%   References: rayi90 wi94

%   AUTHOR: Peter L. Søndergaard
%   TESTING: TEST_PUREFREQ
%   REFERENCE: REF_DSTII

complainif_argnonotinrange(nargin,1,3,mfilename);

if nargin<3
  dim=[];
end;

if nargin<2
  L=[];
end;
   
[f,L,Ls,W,dim,permutedsize,order]=assert_sigreshape_pre(f,L,dim,'DSTII');
 
if ~isempty(L)
  f=postpad(f,L);
end;

c = comp_dst(f,2);

c=assert_sigreshape_post(c,dim,permutedsize,order);

% This is a slow, but convenient way of expressing the above algorithm.
%R=1/sqrt(2)*[zeros(1,L); ...
%	     diag(exp((1:L)*pi*i/(2*L)));...	     
%	     [flipud(diag(-exp(-(1:L-1)*pi*i/(2*L)))),zeros(L-1,1)]];
%R(L+1,L)=i;

%c=i*(R'*fft([f;-flipud(f)])/sqrt(L)/2);


function c=dctii(f,L,dim)
%DCTII  Discrete Consine Transform type II
%   Usage:  c=dctii(f);
%           c=dctii(f,L);
%           c=dctii(f,[],dim);
%           c=dctii(f,L,dim);
%
%   `dctii(f)` computes the discrete cosine transform of type II of the
%   input signal *f*. If *f* is multi-dimensional, the transformation is
%   applied along the first non-singleton dimension.
%
%   `dctii(f,L)` zero-pads or truncates *f* to length *L* before doing the
%   transformation.
%
%   `dctii(f,[],dim)` or `dctii(f,L,dim)` applies the transformation along
%   dimension *dim*.
%
%   The transform is real (output is real if input is real) and orthonormal.
%
%   This is the inverse of |dctiii|.
%
%   Let *f* be a signal of length *L*, let `c=dctii(f)` and define the
%   vector *w* of length *L* by
%
%   .. w = [1/sqrt(2) 1 1 1 1 ...]
%
%   .. math:: w\left(n\right)=\begin{cases}\frac{1}{\sqrt{2}} & \text{if }n=0\\1 & \text{otherwise}\end{cases}
%
%   Then 
%
%   ..                     L-1
%     c(n+1) = sqrt(2/L) * sum w(n+1)*f(m+1)*cos(pi*n*(m+.5)/L) 
%                          m=0 
%
%   .. math:: c\left(n+1\right)=\sqrt{\frac{2}{L}}\sum_{m=0}^{L-1}w\left(n\right)f\left(m+1\right)\cos\left(\frac{\pi}{L} n\left(m+\frac{1}{2}\right)\right)
%
%   Examples:
%   ---------
%
%   The following figures show the first 4 basis functions of the DCTII of
%   length 20:::
%
%     % The dctiii is the adjoint of dctii.
%     F=dctiii(eye(20));
%
%     for ii=1:4
%       subplot(4,1,ii);
%       stem(F(:,ii));
%     end;
%
%   See also:  dctiii, dctiv, dstii
%
%   References: rayi90 wi94

%   AUTHOR: Peter L. Søndergaard
%   TESTING: TEST_PUREFREQ
%   REFERENCE: REF_DCTII

complainif_argnonotinrange(nargin,1,3,mfilename);

if nargin<3
  dim=[];
end;

if nargin<2
  L=[];
end;

[f,L,Ls,W,dim,permutedsize,order]=assert_sigreshape_pre(f,L,dim,'DCTII');

if ~isempty(L)
  f=postpad(f,L);
end;
  c=comp_dct(f,2);
% c=zeros(L,W,assert_classname(f));
% 
% m1=1/sqrt(2)*exp(-(0:L-1)*pi*i/(2*L)).';
% m1(1)=1;
% 
% m2=1/sqrt(2)*exp((1:L-1)*pi*i/(2*L)).';
% 
% s1=fft([f;flipud(f)]);
% 
% % This could be done by a repmat instead.
% for w=1:W
%   c(:,w)=s1(1:L,w).*m1+[0;s1(2*L:-1:L+2,w).*m2];
% end;
% 
% c=c/sqrt(L)/2;

if isreal(f)
  c=real(c);
end;

c=assert_sigreshape_post(c,dim,permutedsize,order);

% This is a slow, but convenient way of expressing the algorithm.
%R=1/sqrt(2)*[diag(exp((0:L-1)*pi*i/(2*L)));...
%	     zeros(1,L); ...
%	     [zeros(L-1,1),flipud(diag(exp(-(1:L-1)*pi*i/(2*L))))]];

%R(1,1)=1;

%c=R'*fft([f;flipud(f)])/sqrt(L)/2;

function c=dctiv(f,L,dim)
%DCTIV  Discrete Consine Transform type IV
%   Usage:  c=dctiv(f);
%
%   `dctiv(f)` computes the discrete cosine transform of type IV of the
%   input signal *f*. If *f* is multi-dimensional, the transformation is
%   applied along the first non-singleton dimension.
%
%   `dctiv(f,L)` zero-pads or truncates *f* to length *L* before doing the
%   transformation.
%
%   `dctiv(f,[],dim)` or `dctiv(f,L,dim)` applies the transformation along
%   dimension *dim*.
%
%   The transform is real (output is real if input is real) and
%   orthonormal.  It is its own inverse.
%
%   Let f be a signal of length *L* and let `c=dctiv(f)`. Then
%
%   ..                     L-1
%     c(n+1) = sqrt(2/L) * sum f(m+1)*cos(pi*(n+.5)*(m+.5)/L) 
%                          m=0 
%
%   .. math:: c\left(n+1\right)=\sqrt{\frac{2}{L}}\sum_{m=0}^{L-1}f\left(m+1\right)\cos\left(\frac{\pi}{L}\left(n+\frac{1}{2}\right)\left(m+\frac{1}{2}\right)\right)
%
%   Examples:
%   ---------
%
%   The following figures show the first 4 basis functions of the DCTIV of
%   length 20:::
%
%     % The dctiv is its own adjoint.
%     F=dctiv(eye(20));
%
%     for ii=1:4
%       subplot(4,1,ii);
%       stem(F(:,ii));
%     end;
%
%   See also:  dctii, dctiii, dstii
%
%   References: rayi90 wi94

%   AUTHOR: Peter L. Søndergaard
%   TESTING: TEST_PUREFREQ
%   REFERENCE: REF_DCTIV

complainif_argnonotinrange(nargin,1,3,mfilename);

if nargin<3
  dim=[];
end;

if nargin<2
  L=[];
end;

[f,L,Ls,W,dim,permutedsize,order]=assert_sigreshape_pre(f,L,dim,'DCTIV');

if ~isempty(L)
  f=postpad(f,L);
end;
c = comp_dct(f,4);
% s1=zeros(2*L,W,assert_classname(f));
% c=zeros(L,W,assert_classname(f));
% 
% m1=1/sqrt(2)*exp(-(0:L-1)*pi*i/(2*L)).';
% m2=1/sqrt(2)*exp((1:L)*pi*i/(2*L)).';
% 
% for w=1:W
%   s1(:,w)=[m1.*f(:,w);flipud(m2).*f(L:-1:1,w)];
% end;
%   
% s1=exp(-pi*i/(4*L))*fft(s1)/sqrt(2*L);
% 
% % This could be done by a repmat instead.
% for w=1:W
%   c(:,w)=s1(1:L,w).*m1+s1(2*L:-1:L+1,w).*m2;
% end;
% 
% if isreal(f)
%   c=real(c);
% end;

c=assert_sigreshape_post(c,dim,permutedsize,order);

% This is a slow, but convenient way of expressing the algorithm.
%R=1/sqrt(2)*[diag(exp(-(0:L-1)*pi*i/(2*L)));...
%	     flipud(diag(exp((1:L)*pi*i/(2*L))))];
  
%c=exp(-pi*i/(4*L))*R.'*fft(R*f)/sqrt(2*L);

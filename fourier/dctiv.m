function c=dctiv(f,L,dim)
%DCTIV  Discrete Consine Transform type IV
%   Usage:  c=dctiv(f);
%
%   DCTIV(f) computes the discrete consine transform of type IV of the
%   input signal f. If f is a matrix, then the transformation is applied to
%   each column.
%
%   DCTIV(f,L) zero-pads or truncates f to length L before doing the
%   transformation.
%
%   DCTIV(f,[],dim) applies the transformation along dimension dim. 
%   DCTIV(f,L,dim) does the same, but pads or truncates to length L.

%   The transform is real (output is real if input is real) and
%   it is orthonormal. It is is own inverse.
%
%   Let f be a signal of length _L and let c=DCTIV(f). Then
%M
%M                         L-1
%M    c(n+1) = sqrt(2/L) * sum f(m+1)*cos(pi*n*(m+.5)/L) 
%M                         m=0 
%F  \[
%F  c\left(n+1\right)=\sqrt{\frac{2}{L}}\sum_{m=0}^{L-1}f\left(m+1\right)\cos\left(\frac{\pi}{L}\left(n+\frac{1}{2}\right)\left(m+\frac{1}{2}\right)\right)
%F  \]
%M
%   See also:  dctii, dctiii, dstii
%M
%R  rayi90 wi94

%   AUTHOR: Peter Soendergaard
%   TESTING: TEST_PUREFREQ
%   REFERENCE: REF_DCTIV

error(nargchk(1,3,nargin));

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

s1=zeros(2*L,W);
c=zeros(L,W);

m1=1/sqrt(2)*exp(-(0:L-1)*pi*i/(2*L)).';
m2=1/sqrt(2)*exp((1:L)*pi*i/(2*L)).';

for w=1:W
  s1(:,w)=[m1.*f(:,w);flipud(m2).*f(L:-1:1,w)];
end;
  
s1=exp(-pi*i/(4*L))*fft(s1)/sqrt(2*L);

% This could be done by a repmat instead.
for w=1:W
  c(:,w)=s1(1:L,w).*m1+s1(2*L:-1:L+1,w).*m2;
end;

if isreal(f)
  c=real(c);
end;

c=assert_sigreshape_post(c,dim,permutedsize,order);

% This is a slow, but convenient way of expressing the algorithm.
%R=1/sqrt(2)*[diag(exp(-(0:L-1)*pi*i/(2*L)));...
%	     flipud(diag(exp((1:L)*pi*i/(2*L))))];
  
%c=exp(-pi*i/(4*L))*R.'*fft(R*f)/sqrt(2*L);

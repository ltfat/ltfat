function c=dsti(f,L,dim)
%DSTI  Discrete Sine Transform type I
%   Usage:  c=dsti(f);
%           c=dsti(f,L);
%           c=dsti(f,[],dim);
%           c=dsti(f,L,dim);
%
%   DSTI(f) computes the discrete sine transform of type I of the
%   input signal f. If f is a matrix, then the transformation is applied to
%   each column. For N-D arrays, the transformation is applied to the first
%   dimension.
%
%   DSTI(f,L) zero-pads or truncates f to length N before doing the
%   transformation.
%
%   DSTI(f,[],dim) applies the transformation along dimension dim. 
%   DSTI(f,L,dim) does the same, but pads or truncates to length L.
%
%   The transform is real (output is real if input is real) and
%   it is orthonormal.
%
%   This transform is its own inverse.
%
%   Let f be a signal of length _L and let c=DSTI(f). Then 
%M
%M                             L-1
%M    c(n+1) = sqrt(2/(L+1)) * sum sin(pi*(n+1)*(m+1)/(L+1)) 
%M                             m=0 
%F  \[
%F  c\left(n+1\right)=\sqrt{\frac{2}{L+1}}\sum_{m=0}^{L-1}f\left(m+1\right)\sin\left(\frac{\pi \left(n+1\right)\left(m+1\right)}{L+1}\right)
%F  \]
%
%   The implementation of this functions uses a simple algorithm that require
%   an FFT of length 2N+2, which might potentially be the product of a large
%   prime number. This may cause the function to sometimes execute slowly.
%   If guaranteed high speed is a concern, please consider using one of the
%   other DST transforms.
%M
%   See also:  dcti, dstiii, dstiv
%M
%R  rayi90 wi94

%   AUTHOR: Peter Soendergaard
%   TESTING: TEST_PUREFREQ
%   REFERENCE: REF_DSTI

  
error(nargchk(1,3,nargin));

if nargin<3
  dim=[];
end;

if nargin<2
  L=[];
end;

[f,L,Ls,W,dim,permutedsize,order]=assert_sigreshape_pre(f,L,dim,'DSTI');

if ~isempty(L)
  f=postpad(f,L);
end;

if L==1
  c=f;
 
else

  c=zeros(L,W);

  s1=dft([zeros(1,W);...
	      f;...
	      zeros(1,W);...
	      -flipud(f)]);


  % This could be done by a repmat instead.
  for w=1:W
    c(:,w)=s1(2:L+1,w)-s1(2*L+2:-1:L+3,w);
  end;

  c=c*i/2;
  
  if isreal(f)
    c=real(c);
  end;

end;

c=assert_sigreshape_post(c,dim,permutedsize,order);

function c=zak(f,a);
%ZAK  Zak transform
%   Usage:  c=zak(f,a);
%
%   ZAK(f,a) computes the Zak transform of f with parameter _a.
%   The coefficients are arranged in an _a x _L/a matrix, where _L is the
%   length of f.
%
%   If f is a matrix, then the transformation is applied to each column.
%   This is then indexed by the third dimension of the output.
%
%
%   Assume that c=ZAK(f,a), where f is a column vector of length L and
%   N=L/a. Then the following holds for $m=0,...,a-1$ and $n=0,...,N-1$
%
%M                         N-1  
%M    c(m+1,n+1)=1/sqrt(N)*sum f(m-k*a+1)*exp(2*pi*i*n*k/N)
%M                         k=0
%F  \begin{eqnarray*}
%F  c(m+1,n+1) & = & \frac{1}{\sqrt{N}}\sum_{k=0}^{N-1}f(m-ka+1)e^{2\pi ink/M}
%F  \end{eqnarray*}
%
%   See also:  izak
%
%R  ja94-4 bohl97-1

%   AUTHOR : Peter Soendergaard
%   TESTING: TEST_ZAK
%   REFERENCE: REF_ZAK

error(nargchk(2,2,nargin));

if (prod(size(a))~=1 || ~isnumeric(a))
  error([callfun,': a must be a scalar']);
end;

if rem(a,1)~=0
  error([callfun,': a must be an integer']);
end;


if size(f,2)>1 && size(f,1)==1
  % f was a row vector.
  f=f(:);
end;

L=size(f,1);
W=size(f,2);
N=L/a;

if rem(N,1)~=0
  error('The parameter for ZAK must divide the length of the signal.');
end;

c=zeros(a,N,W);

for ii=1:W
  % Compute it, it can be done in one line!
  % We use a normalized DFT, as this gives the correct normalization
  % of the Zak transform.
  c(:,:,ii)=dft(reshape(f(:,ii),a,N),[],2);
end;


%OLDFORMAT

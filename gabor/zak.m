function c=zak(f,a);
%ZAK  Zak transform
%   Usage:  c=zak(f,a);
%
%   `zak(f,a)` computes the Zak transform of *f* with parameter *a*.  The
%   coefficients are arranged in an $a \times L/a$ matrix, where *L* is the
%   length of *f*.
%
%   If *f* is a matrix then the transformation is applied to each column.
%   This is then indexed by the third dimension of the output.
%
%   Assume that $c=zak(f,a)$, where *f* is a column vector of length *L* and
%   $N=L/a$. Then the following holds for $m=0,\ldots,a-1$ and $n=0,\ldots,N-1$
%
%   ..                     N-1  
%     c(m+1,n+1)=1/sqrt(N)*sum f(m-k*a+1)*exp(2*pi*i*n*k/N)
%                          k=0
%
%   .. math:: c(m+1,n+1) = \frac{1}{\sqrt{N}}\sum_{k=0}^{N-1}f(m-ka+1)e^{2\pi ink/M}
%
%   Examples:
%   ---------
%
%   This figure shows the absolute value of the Zak-transform of a Gaussian.
%   Notice that the Zak-transform is 0 in only a single point right in the
%   middle of the plot :::
%
%     a=64;
%     L=a^2; 
%     g=pgauss(L);
%     zg=zak(g,a);
%
%     surf(abs(zg));
%   
%   This figure shows the absolute value of the Zak-transform of a 4th order
%   Hermite function.  Notice how the Zak transform of the Hermite functions
%   is zero on a circle centered on the corner :::
%
%     a=64;
%     L=a^2; 
%     g=pherm(L,4);
%     zg=zak(g,a);
%
%     surf(abs(zg));
%
%   See also:  izak
%
%   References: ja94-4 bohl97-1

%   AUTHOR : Peter L. Søndergaard
%   TESTING: TEST_ZAK
%   REFERENCE: REF_ZAK

complainif_argnonotinrange(nargin,2,2,mfilename);

if (prod(size(a))~=1 || ~isnumeric(a))
  error([mfilename,': a must be a scalar']);
end;

if rem(a,1)~=0
  error([mfilename,': a must be an integer']);
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

c=zeros(a,N,W,assert_classname(f));

for ii=1:W
  % Compute it, it can be done in one line!
  % We use a normalized DFT, as this gives the correct normalization
  % of the Zak transform.
  c(:,:,ii)=dft(reshape(f(:,ii),a,N),[],2);
end;



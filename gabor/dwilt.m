function [c,Ls,g]=dwilt(f,g,M,L)
%DWILT  Discrete Wilson transform
%   Usage:  c=dwilt(f,g,M);
%           c=dwilt(f,g,M,L);
%           [c,Ls]=dwilt(...);
%
%   Input parameters:
%         f     : Input data
%         g     : Window function.
%         M     : Number of bands.
%         L     : Length of transform to do.
%   Output parameters:
%         c     : $2M \times N$ array of coefficients.
%         Ls    : Length of input signal.
%
%   `dwilt(f,g,M)` computes a discrete Wilson transform with *M* bands and
%   window *g*.
%
%   The length of the transform will be the smallest possible that is
%   larger than the signal. *f* will be zero-extended to the length of the 
%   transform. If *f* is a matrix, the transformation is applied to each column.
%
%   The window *g* may be a vector of numerical values, a text string or a
%   cell array. See the help of |wilwin|_ for more details.
%
%   `dwilt(f,g,M,L)` computes the Wilson transform as above, but does a
%   transform of length *L*. *f* will be cut or zero-extended to length *L*
%   before the transform is done.
%
%   `[c,Ls]=dwilt(f,g,M)` or `[c,Ls]=dwilt(f,g,M,L)` additionally return the
%   length of the input signal *f*. This is handy for reconstruction::
%
%     [c,Ls]=dwilt(f,g,M);
%     fr=idwilt(c,gd,M,Ls);
%
%   will reconstruct the signal *f* no matter what the length of *f* is, provided
%   that *gd* is a dual Wilson window of *g*.
%
%   `[c,Ls,g]=dwilt(...)` additionally outputs the window used in the
%   transform. This is useful if the window was generated from a description
%   in a string or cell array.
%
%   A Wilson transform is also known as a maximally decimated, even-stacked
%   cosine modulated filter bank.
%
%   Use the function |wil2rect|_ to visualize the coefficients or to work
%   with the coefficients in the TF-plane.
%
%   Assume that the following code has been executed for a column vector *f*::
%
%     c=dwilt(f,g,M);  % Compute a Wilson transform of f.
%     N=size(c,2)*2;   % Number of translation coefficients.
%
%   The following holds for $m=0,\ldots,M-1$ and $n=0,\ldots,N/2-1$:
%
%   If $m=0$:
%
%   ..               L-1 
%     c(m+1,n+1)   = sum f(l+1)*g(l-2*n*M+1)
%                    l=0  
%
%   .. math:: c\left(1,n+1\right) = \sum_{l=0}^{L-1}f(l+1)g\left(l-2nM+1\right)
%
%
%   If $m$ is odd and less than $M$
%
%   ..               L-1 
%     c(m+1,n+1)   = sum f(l+1)*sqrt(2)*sin(pi*m/M*l)*g(k-2*n*M+1)
%                    l=0  
% 
%                    L-1 
%     c(m+M+1,n+1) = sum f(l+1)*sqrt(2)*cos(pi*m/M*l)*g(k-(2*n+1)*M+1)
%                    l=0  
%
%   .. math:: c\left(m+1,n+1\right) & = & \sqrt{2}\sum_{l=0}^{L-1}f(l+1)\sin(\pi\frac{m}{M}l)g(l-2nM+1)\\
%             c\left(m+M+1,n+1\right) & = &
%             \sqrt{2}\sum_{l=0}^{L-1}f(l+1)\cos(\pi\frac{m}{M}l)g\left(l-\left(2n+1\right)M+1\right)
%
%   If $m$ is even and less than $M$
%
%   ..               L-1 
%     c(m+1,n+1)   = sum f(l+1)*sqrt(2)*cos(pi*m/M*l)*g(l-2*n*M+1)
%                    l=0  
% 
%                    L-1 
%     c(m+M+1,n+1) = sum f(l+1)*sqrt(2)*sin(pi*m/M*l)*g(l-(2*n+1)*M+1)
%                    l=0  
%
%   .. math:: c\left(m+1,n+1\right) & = & \sqrt{2}\sum_{l=0}^{L-1}f(l+1)\cos(\pi\frac{m}{M}l)g(l-2nM+1)\\
%             c\left(m+M+1,n+1\right) & = &
%             \sqrt{2}\sum_{l=0}^{L-1}f(l+1)\sin(\pi\frac{m}{M}l)g\left(l-\left(2n+1\right)M+1\right)
%
%   if $m=M$ and $M$ is even:
%
%   ..               L-1 
%     c(m+1,n+1)   = sum f(l+1)*(-1)^(l)*g(l-2*n*M+1)
%                    l=0
%
%   .. math:: c\left(M+1,n+1\right) = \sum_{l=0}^{L-1}f(l+1)(-1)^{l}g(l-2nM+1)
%
%   else if $m=M$ and $M$ is odd
%
%   ..               L-1 
%     c(m+1,n+1)   = sum f(l+1)*(-1)^l*g(l-(2*n+1)*M+1)
%                    l=0
%
%   .. math:: c\left(M+1,n+1\right) = \sum_{k=0}^{L-1}f(l+1)(-1)^{l}g\left(l-\left(2n+1\right)M+1\right)
%
%   See also:  idwilt, wilwin, wil2rect, dgt, wmdct, wilorth
%
%   References: bofegrhl96-1 liva95 dajajo91

%   AUTHOR : Peter L. Søndergaard.
%   TESTING: TEST_DWILT
%   REFERENCE: REF_DWILT

error(nargchk(3,4,nargin));

if nargin<4
  L=[];
end;


assert_squarelat(M,M,1,'DWILT',0);

if ~isempty(L)
  if (prod(size(L))~=1 || ~isnumeric(L))
    error('%s: L must be a scalar','DWILT');
  end;
  
  if rem(L,1)~=0
    error('%s: L must be an integer','DWILT');
  end;
end;

% Change f to correct shape.
[f,Ls,W,wasrow,remembershape]=comp_sigreshape_pre(f,'DWILT',0);

if isempty(L)
  % Smallest length transform.
  Lsmallest=2*M;

  % Choose a transform length larger than the signal
  L=ceil(Ls/Lsmallest)*Lsmallest;
else

  if rem(L,2*M)~=0
    error('%s: The length of the transform must be divisable by 2*M = %i',...
          'DWILT',2*M);
  end;

end;

[g,info]=wilwin(g,M,L,'DWILT');

f=postpad(f,L);

% If the signal is single precision, make the window single precision as
% well to avoid mismatches.
if isa(f,'single')
  g=single(g);
end;

% Call the computational subroutines.
c=comp_dwilt(f,g,M,L);


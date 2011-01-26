function [c,Ls]=wmdct(f,g,M,L)
%WMDCT  Windowed MDCT transform.
%   Usage:  c=wmdct(f,g,M);
%           c=wmdct(f,g,M,L);
%           [c,Ls]=wmdct(f,g,M);
%           [c,Ls]=wmdct(f,g,M,L);
%
%   Input parameters:
%         f     : Input data
%         g     : Window function.
%         M     : Number of bands.
%         L     : Length of transform to do.
%   Output parameters:
%         c     : 2*M x N array of coefficients.
%         Ls    : Length of input signal.
%
%   WMDCT(f,g,M) computes a Windowed Modified Discrete Cosine Transform with
%   M bands and window g.
%
%   The length of the transform will be the smallest possible that is
%   larger than the signal. f will be zero-extended to the length of the 
%   transform. If f is a matrix, the transformation is applied to each column.
%   g must be whole-point even.
%
%   The window g may be a vector of numerical values, a text string or a
%   cell array. See the help of WILWIN for more details.
%
%   WMDCT(f,g,M,L) computes the MDCT transform as above, but does
%   a transform of length L. f will be cut or zero-extended to length L
%   before the transform is done.
%
%   [c,Ls]=WMDCT(f,g,M) or [c,Ls]=WMDCT(f,g,M,L) additionally return the
%   length of the input signal f. This is handy for reconstruction:
%
%C           [c,Ls]=wmdct(f,g,M);
%C           fr=iwmdct(c,gd,M,Ls);
%
%   will reconstuct the signal f no matter what the length of f is, provided
%   that _gd is a dual Wilson window of g.
%
%   The WMDCT is sometimes known as an odd-stacked cosine modulated filter
%   bank. The WMDCT defined by this routine is slightly different from the
%   most common definition of the WMDCT, in order to be able to use the
%   same window functions as the Wilson transform.
%
%   Assume that the following code has been executed for a column vector f
%   of length L:
%
%C           c=wmdct(f,g,M);  % Compute the WMDCT of f.
%C           N=size(c,2);    % Number of translation coefficients.
%
%   The following holds for $m=0,...,M-1$ 
%   and $n=0,...,N-1$:
%
%   If $m+n$ is even:
%
%M                    L-1
%M       c(m+1,n+1) = sum f(l+1)*cos(pi*(m+.5)*l/M+pi/4)*g(l-n*M+1)
%M                    l=0
%F  \begin{eqnarray*}
%F  c\left(m+1,n+1\right) & = & \sqrt{2}\sum_{l=0}^{L-1}f(l+1)\cos\left(\frac{\pi}{M}\left(m+\frac{1}{2}\right)l+\frac{\pi}{4}\right)g(l-nM+1)
%F  \end{eqnarray*}
%
%   If $m+n$ is odd:
%M                    L-1
%M       c(m+1,n+1) = sum f(l+1)*sin(pi*(m+.5)*l/M+pi/4)*g(l-n*M+1)
%M                    l=0
%F  \begin{eqnarray*}
%F  c\left(m+1,n+1\right) & = & \sqrt{2}\sum_{l=0}^{L-1}f(l+1)\sin\left(\frac{\pi}{M}\left(m+\frac{1}{2}\right)l+\frac{\pi}{4}\right)g(l-nM+1)
%F  \end{eqnarray*}
%
%   See also:  iwmdct, wilwin, dwilt, wildual, wilorth
%
%R  prbr86 prjobr87 ma92 bohl96-1 

%   AUTHOR:    Peter Soendergaard
%   TESTING:   TEST_WMDCT
%   REFERENCE: REF_WMDCT

error(nargchk(3,4,nargin));

if nargin<4
  L=[];
end;

[f,g,L,Ls,W,info] = gabpars_from_windowsignal(f,g,M,2*M,L,'WMDCT');

c  = comp_dwiltiii(f,g,M,L);

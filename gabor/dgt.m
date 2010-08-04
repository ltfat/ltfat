function [c,Ls]=dgt(f,g,a,M,varargin)
%DGT  Discrete Gabor transform.
%   Usage:  c=dgt(f,g,a,M);
%           c=dgt(f,g,a,M,L);
%           [c,Ls]=dgt(f,g,a,M);
%           [c,Ls]=dgt(f,g,a,M,L);
%
%   Input parameters:
%         f     : Input data
%         g     : Window function.
%         a     : Length of time shift.
%         M     : Number of channels.
%         L     : Length of transform to do.
%   Output parameters:
%         c     : M*N array of coefficients.
%         Ls    : Length of input signal.
%
%   DGT(f,g,a,M) computes the Gabor coefficients of the input
%   signal f with respect to the window g and parameters _a and M. The
%   output is a vector/matrix in a rectangular layout.
%
%   The length of the transform will be the smallest multiple of _a and M
%   that is larger than the signal. f will be zero-extended to the length of
%   the transform. If f is a matrix, the transformation is applied to each
%   column. The length of the transform done can be obtained by
%   L=size(c,2)*a;
%
%   The window g may be a vector of numerical values, a text string or a
%   cell array. See the help of GABWIN for more detailts.
%
%   DGT(f,g,a,M,L) computes the Gabor coefficients as above, but does
%   a transform of length L. f will be cut or zero-extended to length L before
%   the transform is done.
%
%   [c,Ls]=DGT(f,g,a,M) or [c,Ls]=DGT(f,g,a,M,L) additionally returns the
%   length of the input signal f. This is handy for reconstruction:
%
%C               [c,Ls]=dgt(f,g,a,M);
%C               fr=idgt(c,gd,a,Ls);
%
%   will reconstruct the signal f no matter what the length of f is, provided
%   that _gd is a dual window of g.
%
%   The Discrete Gabor Transform is defined as follows: Consider a window g
%   and a one-dimensional signal f of length L and define $N=L/a$.
%   The output from c=DGT(f,g,a,M) is then given by
%
%M                 L-1 
%M    c(m+1,n+1) = sum f(l+1)*exp(-2*pi*i*m*l/M)*conj(g(l-a*n+1)), 
%M                 l=0  
%F  \[c\left(m+1,n+1\right)=\sum_{l=0}^{L-1}f(l+1)e^{-2\pi ilm/M}\overline{g(l-an+1)}\]
%
%   where $m=0,...,M-1$ and $n=0,...,N-1$ and $l-an$ is computed modulo
%   L.
%
%   DGT takes the following flags at the end of the line of input
%   arguments:
%
%-     'freqinv'  - Compute a DGT using a frequency-invariant phase. This
%                   is the default convention described above.
%
%-     'timeinv'  - Compute a DGT using a time-invariant phase. This
%                   convention is typically used in filter bank algorithms.
%
%   See also:  idgt, gabwin, dwilt, gabdual, phaselock
%
%   Demos:  demo_dgt
%
%R  fest98 gr01

%   AUTHOR : Peter Soendergaard.
%   TESTING: TEST_DGT
%   REFERENCE: REF_DGT
  
% Assert correct input.

if nargin<4
  error('%s: Too few input parameters.',upper(mfilename));
end;

definput.keyvals.L=[];
definput.flags.phase={'freqinv','timeinv'};
[flags,kv]=ltfatarghelper({'L'},definput,varargin);

[f,g,L,Ls] = gabpars_from_windowsignal(f,g,a,M,kv.L);

c=comp_dgt(f,g,a,M,L);

if flags.do_timeinv
  c=phaselock(c,a);
end;


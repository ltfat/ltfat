function [c,Ls,g]=dgt(f,g,a,M,varargin)
%DGT  Discrete Gabor transform
%   Usage:  c=dgt(f,g,a,M);
%           c=dgt(f,g,a,M,L);
%           [c,Ls]=dgt(f,g,a,M);
%           [c,Ls]=dgt(f,g,a,M,L);
%
%   Input parameters:
%         f     : Input data.
%         g     : Window function.
%         a     : Length of time shift.
%         M     : Number of channels.
%         L     : Length of transform to do.
%   Output parameters:
%         c     : $M \times N$ array of coefficients.
%         Ls    : Length of input signal.
%
%   `dgt(f,g,a,M)` computes the Gabor coefficients (also known as a windowed
%   Fourier transform) of the input signal *f* with respect to the window
%   *g* and parameters *a* and *M*. The output is a vector/matrix in a
%   rectangular layout.
%
%   The length of the transform will be the smallest multiple of *a* and *M*
%   that is larger than the signal. *f* will be zero-extended to the length of
%   the transform. If *f* is a matrix, the transformation is applied to each
%   column. The length of the transform done can be obtained by
%   `L=size(c,2)*a;`
%
%   The window *g* may be a vector of numerical values, a text string or a
%   cell array. See the help of |gabwin|_ for more details.
%
%   `dgt(f,g,a,M,L)` computes the Gabor coefficients as above, but does
%   a transform of length *L*. f will be cut or zero-extended to length *L* before
%   the transform is done.
%
%   `[c,Ls]=dgt(f,g,a,M)` or `[c,Ls]=dgt(f,g,a,M,L)` additionally returns the
%   length of the input signal *f*. This is handy for reconstruction::
%
%               [c,Ls]=dgt(f,g,a,M);
%               fr=idgt(c,gd,a,Ls);
%
%   will reconstruct the signal *f* no matter what the length of *f* is, provided
%   that *gd* is a dual window of *g*.
%
%   `[c,Ls,g]=dgt(...)` additionally outputs the window used in the
%   transform. This is useful if the window was generated from a description
%   in a string or cell array.
%
%   The Discrete Gabor Transform is defined as follows: Consider a window *g*
%   and a one-dimensional signal *f* of length *L* and define $N=L/a$.
%   The output from `c=dgt(f,g,a,M)` is then given by:
%
%   ..              L-1 
%      c(m+1,n+1) = sum f(l+1)*conj(g(l-a*n+1))*exp(-2*pi*i*m*l/M), 
%                   l=0  
%
%   .. math:: c\left(m+1,n+1\right)=\sum_{l=0}^{L-1}f(l+1)\overline{g(l-an+1)}e^{-2\pi ilm/M}
%
%   where $m=0,\ldots,M-1$ and $n=0,\ldots,N-1$ and $l-an$ is computed modulo *L*.
%
%   `dgt` takes the following flags at the end of the line of input
%   arguments:
%
%     'freqinv'  Compute a DGT using a frequency-invariant phase. This
%                is the default convention described above.
%
%     'timeinv'  Compute a DGT using a time-invariant phase. This
%                convention is typically used in filter bank algorithms.
%
%   See also:  idgt, gabwin, dwilt, gabdual, phaselock
%
%   Demos:  demo_dgt
% 
%   References: fest98 gr01

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

c=comp_dgt(f,g,a,M,L,flags.do_timeinv);

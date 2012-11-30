function [f,g]=idgt(coef,g,a,varargin)
%IDGT  Inverse discrete Gabor transform
%   Usage:  f=idgt(c,g,a);
%           f=idgt(c,g,a,Ls);
%           f=idgt(c,g,a,Ls,lt);
%
%   Input parameters:
%         c     : Array of coefficients.
%         g     : Window function.
%         a     : Length of time shift.
%         Ls    : Length of signal.
%         lt    : Lattice type (for non-separable lattices)
%   Output parameters:
%         f     : Signal.
%
%   `idgt(c,g,a)` computes the Gabor expansion of the input coefficients
%   *c* with respect to the window *g* and time shift *a*. The number of 
%   channels is deduced from the size of the coefficients *c*.
%
%   `idgt(c,g,a,Ls)` does as above but cuts or extends *f* to length *Ls*.
%
%   `[f,g]=idgt(...)` additionally outputs the window used in the
%   transform. This is useful if the window was generated from a description
%   in a string or cell array.
%
%   For perfect reconstruction, the window used must be a dual window of the
%   one used to generate the coefficients.
%
%   The window *g* may be a vector of numerical values, a text string or a
%   cell array. See the help of |gabwin|_ for more details.
%
%   If *g* is a row vector, then the output will also be a row vector. If *c* is
%   3-dimensional, then `idgt` will return a matrix consisting of one column
%   vector for each of the TF-planes in *c*.
%
%   Assume that `f=idgt(c,g,a,L)` for an array *c* of size $M \times N$. 
%   Then the following holds for $k=0,\ldots,L-1$: 
% 
%   ..          N-1 M-1          
%     f(l+1)  = sum sum c(m+1,n+1)*exp(2*pi*i*m*l/M)*g(l-a*n+1)
%               n=0 m=0          
%
%   .. math:: f(l+1) = \sum_{n=0}^{N-1}\sum_{m=0}^{M-1}c(m+1,n+1)e^{2\pi
%      iml/M}g(l-an+1)
%
%   Non-separable lattices:
%   -----------------------
%
%   `idgt(c,g,a,'lt',lt)` computes the Gabor expansion of the input
%   coefficients *c* with respect to the window *g*, time shift *a* and
%   lattice type *lt*. Please see the help of |matrix2latticetype|_ for a
%   precise description of the parameter *lt*.
%
%   Assume that `f=dgt(c,g,a,L,lt)` for an array *c* of size $M\times N$.
%   Then the following holds for $k=0,\ldots,L-1$:
% 
%   ..          N-1 M-1          
%     f(l+1)  = sum sum c(m+1,n+1)*exp(2*pi*i*m*l/M)*g(l-a*n+1)
%               n=0 m=0          
%
%   .. math:: f(l+1) = \sum_{n=0}^{N-1}\sum_{m=0}^{M-1}c(m+1,n+1)e^{2\pi iml/M}g(l-an+1)
%
%   Additional parameters:
%   ----------------------
%
%   `idgt` takes the following flags at the end of the line of input
%   arguments:
%
%     'freqinv'  Compute an IDGT using a frequency-invariant phase. This
%                is the default convention described above.
%
%     'timeinv'  Compute an IDGT using a time-invariant phase. This
%                convention is typically used in filter bank algorithms.
%
%   See also:  dgt, gabwin, dwilt, gabtight

%   AUTHOR : Peter L. Søndergaard.
%   TESTING: TEST_DGT
%   REFERENCE: OK

% Check input paramameters.

if nargin<3
  error('%s: Too few input parameters.',upper(mfilename));
end;

if numel(g)==1
  error('g must be a vector (you probably forgot to supply the window function as input parameter.)');
end;

definput.keyvals.Ls=[];
definput.keyvals.lt=[0 1];
definput.flags.phase={'freqinv','timeinv'};
[flags,kv,Ls]=ltfatarghelper({'Ls'},definput,varargin);

M=size(coef,1);
N=size(coef,2);
W=size(coef,3);

if ~isnumeric(a) || ~isscalar(a)
  error('%s: "a" must be a scalar',upper(mfilename));
end;

if rem(a,1)~=0
  error('%s: "a" must be an integer',upper(mfilename));
end;

L=N*a;

Ltest=dgtlength(L,a,M,kv.lt);

if Ltest~=L
    error(['%s: Incorrect size of coefficient array or "a" parameter. See ' ...
           'the help of DGTLENGTH for the requirements.'], ...
          upper(mfilename))
end;

g=gabwin(g,a,M,L,kv.lt,'callfun',upper(mfilename));

f=comp_idgt(coef,g,a,kv.lt,flags.do_timeinv,0);

% Cut or extend f to the correct length, if desired.
if ~isempty(Ls)
  f=postpad(f,Ls);
else
  Ls=L;
end;

f=comp_sigreshape_post(f,Ls,0,[0; W]);

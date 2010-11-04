function [f]=idgt(coef,g,a,varargin)
%IDGT  Inverse discrete Gabor transform.
%   Usage:  f=idgt(c,g,a);
%           f=idgt(c,g,a,Ls);
%
%   Input parameters:
%         c     : Array of coefficients.
%         g     : Window function.
%         a     : Length of time shift.
%         Ls    : length of signal.
%   Output parameters:
%         f     : Signal.
%
%   IDGT(c,g,a) computes the Gabor expansion of the input coefficients
%   c with respect to the window g and time shift _a. The number of 
%   channels is deduced from the size of the coefficients c.
%
%   IDGT(c,g,a,Ls) does as above but cuts or extends f to length Ls.
%
%   For perfect reconstruction, the window used must be a dual window of the
%   one used to generate the coefficients.
%
%   The window g may be a vector of numerical values, a text string or a
%   cell array. See the help of GABWIN for more detailts.
%
%   If g is a row vector, then the output will also be a row vector. If c is
%   3-dimensional, then IDGT will return a matrix consisting of one column
%   vector for each of the TF-planes in c.
%
%   Assume that f=IDGT(c,g,a,L) for an array c of size $M$ x $N$. 
%   Then the following holds for k=0,...,L-1: 
% 
%M            N-1 M-1          
%M  f(l+1)  = sum sum c(m+1,n+1)*exp(2*pi*i*m*l/M)*g(l-a*n+1)
%M            n=0 m=0          
%F  \begin{eqnarray*}
%F  f(l+1) & = & \sum_{n=0}^{N-1}\sum_{m=0}^{M-1}c(m+1,n+1)e^{2\pi iml/M}g(l-an+1)
%F  \end{eqnarray*}
%
%   IDGT takes the following flags at the end of the line of input
%   arguments:
%
%-     'freqinv'  - Compute an IDGT using a frequency-invariant phase. This
%                   is the default convention described above.
%
%-     'timeinv'  - Compute an IDGT using a time-invariant phase. This
%                   convention is typically used in filter bank algorithms.
%
%   See also:  dgt, gabwin, dwilt, gabtight

%   AUTHOR : Peter Soendergaard.
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
definput.flags.phase={'freqinv','timeinv'};
[flags,kv,Ls]=ltfatarghelper({'Ls'},definput,varargin);

wasrow=0;

if isnumeric(g)
  if size(g,2)>1
    if size(g,1)>1
      error('g must be a vector');
    else
      % g was a row vector.
      g=g(:);
      
      % If the input window is a row vector, and the dimension of c is
      % equal to two, the output signal will also
      % be a row vector.
      if ndims(coef)==2
        wasrow=1;
      end;
    end;
  end;
end;

M=size(coef,1);
N=size(coef,2);
W=size(coef,3);

% use assert_squarelat to check a and the window size.
assert_squarelat(a,M,1,'IDGT');

L=N*a;

g=comp_window(g,a,M,L,0,'IDGT');

assert_L(L,size(g,1),L,a,M,'IDGT');

f=comp_idgt(coef,g,a,M,L,flags.do_timeinv);

% Cut or extend f to the correct length, if desired.
if ~isempty(Ls)
  f=postpad(f,Ls);
else
  Ls=L;
end;

f=comp_sigreshape_post(f,Ls,wasrow,[0; W]);





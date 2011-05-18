function [f,g]=idgtreal(coef,g,a,M,varargin)
%IDGTREAL  Inverse discrete Gabor transform.
%   Usage:  f=idgtreal(c,g,a,M);
%           f=idgtreal(c,g,a,M,Ls);
%
%   Input parameters:
%         c     : Array of coefficients.
%         g     : Window function.
%         a     : Length of time shift.
%         M     : Number of channels.
%         Ls    : length of signal.
%   Output parameters:
%         f     : Signal.
%
%   IDGTREAL(c,g,a,M) computes the Gabor expansion of the input coefficients
%   c with respect to the real-valued window g, time shift _a and number of
%   channels M. c is assumed to be the positive frequencies of the Gabor
%   expansion of a real-valued signal.
%
%   It must hold that size(c,1)==floor(M/2)+1. Note that since the
%   correct number of channels cannot be deduced from the input, IDGTREAL
%   takes an additional parameter as opposed to IDGT.
%
%   The window g may be a vector of numerical values, a text string or a
%   cell array. See the help of GABWIN for more details.
%  
%   IDGTREAL(c,g,a,M,Ls) does as above but cuts or extends f to length Ls.
%
%   [f,g]=IDGTREAL(...) additionally outputs the window used in the
%   transform. This is usefull if the window was generated from a description
%   in a string or cell array.
%
%   For perfect reconstruction, the window used must be a dual window of the
%   one used to generate the coefficients.
%
%   If g is a row vector, then the output will also be a row vector. If c is
%   3-dimensional, then IDGTREAL will return a matrix consisting of one column
%   vector for each of the TF-planes in c.
%
%   See the help on IDGT for the precise definition of the inverse Gabor
%   transform.
%
%   IDGTREAL takes the following flags at the end of the line of input
%   arguments:
%
%-     'freqinv'  - Compute an IDGT using a frequency-invariant phase. This
%                   is the default convention described in the help for IDGT.
%
%-     'timeinv'  - Compute an IDGT using a time-invariant phase. This
%                   convention is typically used in filter bank algorithms.
%
%   See also:  idgt, gabwin, gabdual, dwilt

%   AUTHOR : Peter Soendergaard.
%   TESTING: TEST_DGT
%   REFERENCE: OK

% Check input paramameters.

if nargin<4
  error('%s: Too few input parameters.',upper(mfilename));
end;

if prod(size(g))==1
  error('g must be a vector (you probably forgot to supply the window function as input parameter.)');
end;

% Define initial value for flags and key/value pairs.
definput.keyvals.Ls=[];
definput.keyvals.dim=[];
definput.flags.phase={'freqinv','timeinv'};

[flags,kv]=ltfatarghelper({'Ls','dim',},definput,varargin);

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

N=size(coef,2);
W=size(coef,3);
M2=floor(M/2)+1;

if M2~=size(coef,1)
  error('Mismatch between the specified number of channels and the size of the input coefficients.');
end;

% use assert_squarelat to check a and the window size.
assert_squarelat(a,M,1,'IDGTREAL');

L=N*a;

g=gabwin(g,a,M,L,'IDGTREAL');

if ~isreal(g)
  error('Window must be real-valued.');
end;

assert_L(L,size(g,1),L,a,M,'IDGTREAL');

% Do the actual computation.
f=comp_idgtreal(coef,g,a,M,L,flags.do_timeinv);

% Cut or extend f to the correct length, if desired.
if ~isempty(kv.Ls)
  f=postpad(f,kv.Ls);
else
  kv.Ls=L;
end;

f=comp_sigreshape_post(f,kv.Ls,wasrow,[0; W]);





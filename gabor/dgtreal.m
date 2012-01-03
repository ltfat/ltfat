function [c,Ls,g]=dgtreal(f,g,a,M,varargin)
%DGTREAL  Discrete Gabor transform.
%   Usage:  c=dgtreal(f,g,a,M);
%           c=dgtreal(f,g,a,M,L);
%           [c,Ls]=dgtreal(f,g,a,M);
%           [c,Ls]=dgtreal(f,g,a,M,L);
%
%   Input parameters:
%         f     : Input data
%         g     : Window function.
%         a     : Length of time shift.
%         M     : Number of modulations.
%         L     : Length of transform to do.
%   Output parameters:
%         c     : M*N array of coefficients.
%         Ls    : Length of input signal.
%
%   DGTREAL(f,g,a,M) computes the Gabor coefficients of the real-valued
%   input signal f with respect to the real-valued window g and parameters
%   _a and M. The output is a vector/matrix in a rectangular layout.
%
%   As opposed to DGT only the coefficients of the positive frequencies of
%   the output are returned. DGTREAL will refuse to work for complex
%   valued input signals.
%
%   The length of the transform will be the smallest multiple of a and M
%   that is larger than the signal. f will be zero-extended to the length of
%   the transform. If f is a matrix, the transformation is applied to each
%   column. The length of the transform done can be obtained by
%   L=size(c,2)*a;
%
%   The window g may be a vector of numerical values, a text string or a
%   cell array. See the help of GABWIN for more details.
%
%   DGTREAL(f,g,a,M,L) computes the Gabor coefficients as above, but does
%   a transform of length L. f will be cut or zero-extended to length L before
%   the transform is done.
%
%   [c,Ls]=DGTREAL(f,g,a,M) or [c,Ls]=DGTREAL(f,g,a,M,L) additionally
%   returns the length of the input signal f. This is handy for
%   reconstruction:
%
%C               [c,Ls]=dgtreal(f,g,a,M);
%C               fr=idgtreal(c,gd,a,M,Ls);
%
%   will reconstruct the signal f no matter what the length of f is, provided
%   that _gd is a dual window of g.
%
%   [c,Ls,g]=DGTREAL(...) additionally outputs the window used in the
%   transform. This is useful if the window was generated from a description
%   in a string or cell array.
%
%   See the help on DGT for the definition of the discrete Gabor
%   transform. This routine will return the coefficients for channel
%   frequencies from 0 to floor(M/2).
%
%   DGTREAL takes the following flags at the end of the line of input
%   arguments:
%
%-     'freqinv'  - Compute a DGTREAL using a frequency-invariant phase. This
%                   is the default convention described in the help for DGT.
%
%-     'timeinv'  - Compute a DGTREAL using a time-invariant phase. This
%                   convention is typically used in filter bank algorithms.

%
%   See also:  dgt, idgtreal, gabwin, dwilt, gabtight
%
%R  fest98 gr01

%   AUTHOR : Peter Soendergaard.
%   TESTING: TEST_DGT
%   REFERENCE: OK
  
% Assert correct input.

if nargin<4
  error('%s: Too few input parameters.',upper(mfilename));
end;

definput.keyvals.L=[];
definput.flags.phase={'freqinv','timeinv'};
[flags,kv]=ltfatarghelper({'L'},definput,varargin);

[f,g,L,Ls] = gabpars_from_windowsignal(f,g,a,M,kv.L);

if ~isreal(g)
  error('The window must be real-valued.');  
end;

c=comp_dgtreal(f,g,a,M,L,flags.do_timeinv);



%OLDFORMAT

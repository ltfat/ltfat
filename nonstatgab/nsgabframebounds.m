function [AF,BF]=nsgabframebounds(g,a,M)
%NSGABFRAMEBOUNDS  Frame bounds of non-stationary Gabor frame
%   Usage:  fcond=nsgabframebounds(g,a,M);
%           [A,B]=nsgabframebounds(g,a,M);
%
%   Input parameters:
%         g     : Cell array of windows
%         a     : Vector of time positions of windows.
%         M     : Vector of numbers of frequency channels.
%   Output parameters:
%         fcond : Frame condition number ($B/A$)
%         A,B   : Frame bounds.
%
%   `nsgabframebounds(g,a,Ls)` calculates the ratio $B/A$ of the frame
%   bounds of the non-stationary discrete Gabor frame defined by windows
%   given in *g* at positions given by *a*. Please see the help on |nsdgt|
%   for a more thourough description of *g* and *a*.
%
%   `[A,B]=nsgabframebounds(g,a,Ls)` returns the actual frame bounds *A*
%   and *B* instead of just the their ratio.
%
%   The computed frame bounds are only valid for the 'painless case' when
%   the number of frequency channels used for computation of |nsdgt| is greater
%   than or equal to the window length. This correspond to cases for which
%   the frame operator is diagonal.
%
%   See also:  nsgabtight, nsdgt, insdgt
%
%   References: ltfatnote018
  
%   AUTHOR : Florent Jaillet
%   TESTING: TEST_NSDGT

% Compute the diagonal of the frame operator.
f=nsgabframediag(g,a,M);

AF=min(f);
BF=max(f);

if nargout<2
  % Avoid the potential warning about division by zero.
  if AF==0
    AF=Inf;
  else
    AF=BF/AF;
  end;
end;


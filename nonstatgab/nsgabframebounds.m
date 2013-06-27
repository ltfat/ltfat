function [AF,BF]=nsgabframebounds(g,a,Ls)
%NSGABFRAMEBOUNDS  Frame bounds of nonstationnary Gabor frame
%   Usage:  fcond=nsgabframebounds(g,a,Ls);
%           [A,B]=nsgabframebounds(g,a,Ls);
%
%   Input parameters:
%         g     : Cell array of windows
%         a     : Vector of time positions of windows.
%         Ls    : Length of analyzed signal.
%   Output parameters:
%         fcond : Frame condition number ($B/A$)
%         A,B   : Frame bounds.
%
%   `nsgabframebounds(g,a,Ls)` calculates the ratio $B/A$ of the frame
%   bounds of the nonstationary discrete Gabor frame defined by windows
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

timepos=cumsum(a)-a(1);

N=length(a); % Number of time positions

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


function [AF,BF]=nonsepgabframebounds(g,a,M,lt,L)
%NONSEPGABFRAMEBOUNDS  Calculate frame bounds of Gabor frame
%   Usage:  fcond=nonsepgabframebounds(g,a,M,lt);
%           [A,B]=nonsepgabframebounds(g,a,M,lt);
%           [A,B]=nonsepgabframebounds(g,a,M,lt,L);
%
%   Input parameters:
%           g     : The window function.
%           a     : Length of time shift.
%           M     : Number of channels.
%           lt    : Lattice type
%           L     : Length of transform to consider.
%   Output parameters:
%           fcond : Frame condition number (B/A)
%           A,B   : Frame bounds.
%          
%   `nonsepgabframebounds(g,a,M,lt)` calculates the ratio $B/A$ of the frame
%   bounds of the non-separable Gabor system with window *g*, and parameters
%   *a*, *M* and *lt*.
%
%   `[A,B]=nonsepgabframebounds(g,a,M,lt)` returns the frame bounds *A* and
%   *B* instead of just the ratio.
%
%   The window *g* may be a vector of numerical values, a text string or a
%   cell array. See the help of |gabwin|_ for more details.
%  
%   If the optional parameter *L* is specified, the window is cut or
%   zero-extended to length *L*.
%
%   See also: gabrieszbounds, gabwin

error(nargchk(4,5,nargin));

if nargin<5
  L=[];
end;

[g,L] = nonsepgabpars_from_window(g,a,M,lt,L);

g=fir2long(g,L);

% ---- algorithm starts here: Use the multi-window factorization ---

% Convert to multi-window
mwin=comp_nonsepwin2multi(g,a,M,lt);

% Get the factorization of the window.
gf=comp_wfac(mwin,a*lt(2),M);

% Compute all eigenvalues.
lambdas=comp_gfeigs(gf,L,a*lt(2),M);
s=size(lambdas,1);

% Min and max eigenvalue.
if a>M
  % This cannot be a frame, so A is identically 0.
  AF=0;
else
  AF=lambdas(1);
end;

BF=lambdas(s);

if nargout<2
  % Avoid the potential warning about division by zero.
  if AF==0
    AF=Inf;
  else
    AF=BF/AF;
  end;
end;


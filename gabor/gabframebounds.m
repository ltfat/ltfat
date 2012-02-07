function [AF,BF]=gabframebounds(g,a,M,L)
%GABFRAMEBOUNDS  Calculate frame bounds of Gabor frame.
%   Usage:  fcond=gabframebounds(g,a,M);
%           [A,B]=gabframebounds(g,a,M);
%           [A,B]=gabframebounds(g,a,M,L);
%
%   Input parameters:
%           g     : The window function.
%           a     : Length of time shift.
%           M     : Number of channels.
%           L     : Length of transform to consider.
%   Output parameters:
%           fcond : Frame condition number (B/A)
%           A,B   : Frame bounds.
%          
%   `gabframebounds(g,a,M)` calculates the ratio $B/A$ of the frame bounds
%   of the Gabor system with window *g*, and parameters *a*, *M*.
%
%   `[A,B]=gabframebounds(g,a,M)` returns the frame bounds *A* and *B*
%   instead of just the ratio.
%
%   The window *g* may be a vector of numerical values, a text string or a
%   cell array. See the help of |gabwin|_ for more details.
%  
%   If the optional parameter *L* is specified, the window is cut or
%   zero-extended to length *L*.
%
%   See also: gabrieszbounds, gabwin

error(nargchk(3,4,nargin));

if nargin<4
  L=[];
end;

[g,L,info] = gabpars_from_window(g,a,M,L);

g=fir2long(g,L);

% Get the factorization of the window.
gf=comp_wfac(g,a,M);

% Compute all eigenvalues.
lambdas=comp_gfeigs(gf,L,a,M);
s=size(lambdas,1);

% Min and max eigenvalue.
if a>M
  % This can is not a frame, so A is identically 0.
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


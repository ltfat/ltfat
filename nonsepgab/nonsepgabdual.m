function [gd,gdfull,gdmatch]=nonsepgabdual(g,a,M,lt,L)
%NONSEPGABDUAL  Canonical dual window of Gabor frame
%   Usage:  gd=nonsepgabdual(g,a,M,lt);
%           gd=nonsepgabdual(g,a,M,lt,L);
%
%   Input parameters:
%         g     : Gabor window.
%         a     : Length of time shift.
%         M     : Number of channels.
%         lt    : Lattice type
%         L     : Length of window. (optional)
%   Output parameters:
%         gd : Canonical dual window.
%
%   `nonsepgabdual(g,a,M,lt)` computes the canonical dual window of the
%   discrete Gabor frame with window *g* and parameters *a*, *M* and *lt*.
%
%   The window *g* may be a vector of numerical values, a text string or a
%   cell array. See the help of |gabwin|_ for more details.
%
%   If the length of *g* is equal to *M*, then the input window is assumed
%   to be an FIR window. In this case, the canonical dual window also has
%   length of *M*. Otherwise the smallest possible transform length is chosen
%   as the window length.
%
%   `nonsepgabdual(g,a,M,lt,L)` returns a window that is the dual window for a
%   system of length *L*. Unless the dual window is a FIR window, the dual
%   window will have length *L*.
%
%   If $a>M$ then the dual window of the Gabor Riesz sequence with window
%   *g* and parameters *a* and *M* will be calculated.

%   AUTHOR : Peter Soendergaard.
%   TESTING: TEST_DGT
%   REFERENCE: REF_NONSEPGABDUAL.
  
% ------ Direct checks on input parameters ----

error(nargchk(4,5,nargin));  

if nargin==4
  L=[];
end;

[g,L,info] = nonsepgabpars_from_window(g,a,M,lt,L);

% -------- Are we in the Riesz sequence of in the frame case

scale=1;
if a>M
  % Handle the Riesz basis (dual lattice) case.
  % Swap a and M, and scale differently.
  scale=a/M;
  tmp=a;
  a=M;
  M=tmp;
end;

% -------- Compute ------------- 

% Just in case, otherwise the call is harmless. 
g=fir2long(g,L);

mwin=comp_nonsepwin2multi(g,a,M,lt);

gdfull=comp_gabdual_long(mwin,a*lt(2),M)*scale;

% We need just the first vector
gd=gdfull(:,1);

gdmatch=comp_nonsepwin2multi(gd,a,M,lt);

% --------- post process result -------
      
if info.wasrow
  gd=gd.';
end;

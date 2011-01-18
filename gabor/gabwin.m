function [g,info] = gabwin(g,a,M,L,callfun);
%GABWIN  Compute a Gabor window from text or cell array.
%   Usage: [g,info] = gabwin(g,a,M,L);
%
%   [g,info]=GABWIN(g,a,M,L) computes a window that fits well with the
%   specified number of channels M, time shift _a and transform length
%   L. The window itself is specified by a text description or a cell array
%   containing additional parameters.
%
%   The window can be specified directly as a vector of numerical
%   values. In this case, GABWIN only checks assumptions about transform
%   sizes etc.
%
%   [g,info]=GABWIN(g,a,M) does the same, but the window must be a FIR
%   window, as the transform length is unspecified.
%
%   The window can be specified as one of the following text strings:
%  
%-      'gauss'     - Gaussian window fitted to the lattice, i.e. tfr=a*M/L.
%
%-      'dualgauss' - Canonical dual of Gaussian window.
%
%-      'tight'     - Tight window generated from a Gaussian.
%
%   In these cases, a long window is generated with a length of L.
%
%   It is also possible to specify one of the window names from FIRWIN. In
%   such a case, GABWIN will generate the speficied FIR window with a length
%   of M.
%
%   The window can also be specified as cell array. The possibilities are:
%
%-    {'gauss',...} - Additional parameters are passed to PGAUSS
%
%-    {'dual',...}  - Canonical dual window of whatever follows. See the
%                   examples below.
%
%-    {'tight',...} - Canonical tight window of whatever follows.
%
%   It is also possible to specify one of the window names from FIRWIN as
%   the first field in the cell array. In this case, the remaining
%   entries of the cell array are passed directly to FIRWIN.
%
%   Some examples:
%
%C     g=gabwin('gauss',a,M,L);
%
%   This computes a Gaussian window of length L fitted for a system with
%   time-shift _a and M channels.
%
%C    g=gabwin({'gauss',1},a,M,L);
%
%   This computes a Gaussian window with equal time and frequency support
%   irrespective of _a and M.
%
%C    gd=gabwin('gaussdual',a,M,L);
%
%   This computes the canonical dual of a Gaussian window fitted for a
%   system with time-shift _a and M channels.
%
%C    gd=gabwin({'tight','gauss'},a,M,L);
%
%   This computes the canonical tight window of the Gaussian window fitted
%   for the system.
%
%C    g=gabwin({'dual',{'hann',20}},a,M,L);
%
%   This computes the dual of a Hann window of length 20.  
%
%   The structure info provides some information about the computed
%   window:
%
%-     info.gauss   - True if the window is a Gaussian.
%
%-     info.tfr     - Time/frequency support ratio of the window. Set
%                     whenever it makes sense.
%
%-     info.wasrow  - Input was a row window
%
%-     info.isfir   - Input is a FIR window
%
%-     info.isdual  - Output is the dual window of the aux window.
%
%-     info.istight - Output is known to be a tight window.
%
%-     info.auxinfo - Info about auxiliary window.
%
%-     info.gl      - Length of window.
%
%   See also: pgauss, firwin, wilwin
  
% Assert correct input.
error(nargchk(3,5,nargin));

if nargin==3
  L=[];
end;

[g,info] = comp_window(g,a,M,L,'GABWIN');

if (info.isfir)  
  if info.istight
    g=g/sqrt(2);
  end;  
end;
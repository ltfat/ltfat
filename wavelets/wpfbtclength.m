function [Lc,L]=wpfbtclength(Ls,wt,varargin)
%WPFBTCLENGTH  WPFBT subband length from a signal length
%   Usage: Lc=wpfbtclength(Ls,wt);
%          [Lc,L]=wpfbtclength(Ls,wt);
%
%   `Lc=wpfbtclength(Ls,wt)` returns the lengths of coefficient subbands 
%   obtained from |wpfbt| for a signal of length *Ls*. Please see the help 
%   on |wpfbt| for an explanation of the parameter *wt*. 
%
%   `[Lc,L]=wpfbtclength(...)` additionally returns the next legal length 
%   of the input signal for the given extension type.
%
%   The function support the same boundary-handling flags as the |fwt|
%   does.
%
%   See also: wpfbt

% AUTHOR: Zdenek Prusa


definput.import = {'fwt'};
[flags,kv]=ltfatarghelper({},definput,varargin);

% Initialize the wavelet filters structure
wt = wfbtinit(wt);

if(flags.do_per)
   a = treeSub(wt);
   L = filterbanklength(Ls,a);
else
   L = Ls;
end

wtPath = nodeBForder(0,wt);
Lc = nodesOutLen(wtPath,L,[],flags.do_per,wt);



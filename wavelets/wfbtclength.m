function [Lc,L]=wfbtclength(Ls,wt,varargin)
%WFBTCLENGTH  WFBT subband lengths from a signal length
%   Usage: Lc=wfbtclength(Ls,wt);
%          [Lc,L]=wfbtclength(...);
%
%   `Lc=wfbtclength(Ls,wt)` returns the lengths of coefficient subbands 
%   obtained from |wfbt| for a signal of length *Ls*. Please see the help 
%   on |wfbt| for an explanation of the parameters *wt*. 
%
%   `[Lc,L]=wfbtclength(...)` additionally returns the next legal length 
%   of the input signal for the given extension type.
%
%   The function support the same boundary-handling flags as the |fwt|
%   does.
%
%   See also: wfbt, wfbtlength

% AUTHOR: Zdenek Prusa

complainif_notposint(Ls,'Ls','WFBTCLENGTH');

definput.import = {'fwt'};
[flags,kv]=ltfatarghelper({},definput,varargin);

% Initialize the wavelet filters structure
wt = wfbtinit(wt);


if(flags.do_per)
   a = treeSub(wt);
   L = filterbanklength(Ls,a);
   Lc = L./a;
else
   L = Ls;
   Lc = treeOutLen(L,0,wt);
end



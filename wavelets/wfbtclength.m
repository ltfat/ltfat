function [Lc,L]=wfbtclength(Ls,wt,varargin);
%WFBTLENGTH  WFBT subband lengthf from a signal length
%   Usage: L=wfbtlength(Ls,wt);
%          L=wfbtlength(Ls,wt,...);
%
%   `wfbtlength(Ls,wt)` returns the length of a Wavelet system that is long
%   enough to expand a signal of length *Ls*. Please see the help on
%   |wfbt| for an explanation of the parameter *wt*.
%
%   If the returned length is longer than the signal length, the signal
%   will be zero-padded by |wfbt| to length *L*.
%
%   See also: wfbt, fwt


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



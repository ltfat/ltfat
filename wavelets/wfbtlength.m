function L=wfbtlength(Ls,wt,varargin);
%WFBTLENGTH  WFBT length from signal
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

complainif_notposint(Ls,'Ls','WFBTLENGTH');

definput.import = {'fwt'};
[flags,kv]=ltfatarghelper({},definput,varargin);

% Initialize the wavelet filters structure
if ~isstruct(wt)
   wt = wfbtinit(wt);
end

if(flags.do_per)
   a = treeSub(wt);
   L = filterbanklength(Ls,a);
else
   L = Ls;
end

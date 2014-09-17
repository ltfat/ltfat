function L=wfbtlength(Ls,wt,varargin);
%WFBTLENGTH  WFBT length from signal
%   Usage: L=wfbtlength(Ls,wt);
%
%   `wfbtlength(Ls,wt)` returns the length of a Wavelet system that is long
%   enough to expand a signal of length *Ls*. Please see the help on
%   |wfbt| for an explanation of the parameter *wt*.
%
%   If the returned length is longer than the signal length, the signal
%   will be zero-padded by |wfbt| to length *L*.
%
%   In addition, the function accepts flags defining boundary extension
%   technique as in |wfbt|. The returned length can be longer than the
%   signal length only in case of `'per'` (periodic extension).
%
%   See also: wfbt, fwt

% AUTHOR: Zdenek Prusa

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

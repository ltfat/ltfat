function [Lc,L]=wpfbtclength(Ls,wt,varargin)
%WPFBTCLENGTH  WPFBT subband length from a signal length
%   Usage: L=wpfbtclength(Ls,wt);
%          L=wpfbtclength(Ls,wt,...);
%
%   `[Lc,L]=wpfbtclength(Ls,wt)` returns the length *L* of a wavelet system
%   that is long enough to expand a signal of length *Ls* and associated
%   vector subband lengths *Lc*. Please see the help on
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
else
   L = Ls;
end

wtPath = nodesBForder(wt);
Lc = nodeOutLen(wtPath,L,[],flags.do_per,wt);



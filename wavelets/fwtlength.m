function L=fwtlength(Ls,w,J,varargin)
%FWTLENGTH  FWT length from signal
%   Usage: L=fwtlength(Ls,w,J);
%
%   `fwtlength(Ls,w,J)` returns the length of a Wavelet system that is long
%   enough to expand a signal of length *Ls*. Please see the help on
%   |fwt| for an explanation of the parameters *w* and *J*.
%
%   If the returned length is longer than the signal length, the signal
%   will be zero-padded by |fwt| to length *L*.
%
%   In addition, the function accepts flags defining boundary extension
%   technique as in |fwt|. The returned length can be longer than the
%   signal length only in case of `'per'` (periodic extension).
%
%   See also: fwt

% AUTHOR: Zdenek Prusa

complainif_notposint(Ls,'Ls','FWTLENGTH');
complainif_notposint(J,'J','FWTLENGTH');

% Initialize the wavelet filters structure
w = fwtinit(w);

definput.import = {'fwtext'};
[flags,kv]=ltfatarghelper({},definput,varargin);

if flags.do_per
   blocksize=w.a(1)^J;
   L=ceil(Ls/blocksize)*blocksize;
% elseif flags.do_valid
%    m = numel(w.g{1}.h);
%    a = w.a(1);
%    rred = (a^J-1)/(a-1)*(m-a);
%    blocksize=w.a(1)^J;
%    L=rred+floor((Ls-rred)/blocksize)*blocksize;
else
   L = Ls;
end

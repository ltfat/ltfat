function L=fwtlength(Ls,w,J,varargin);
%FWTLENGTH  FWT length from signal
%   Usage: L=fwtlength(Ls,w,J);
%          L=fwtlength(Ls,w,J,...);
%
%   `fwtlength(Ls,w,J)` returns the length of a Wavelet system that is long
%   enough to expand a signal of length *Ls*. Please see the help on
%   |fwt| for an explanation of the parameters *h* and *J*.
%
%   If the returned length is longer than the signal length, the signal
%   will be zero-padded by |fwt| to length *L*.
%
%   See also: fwt

% Initialize the wavelet filters structure
w = fwtinit(w);

definput.import = {'fwtext'};
[flags,kv]=ltfatarghelper({},definput,varargin);

if flags.do_per
   blocksize=w.a(1)^J;
   L=ceil(Ls/blocksize)*blocksize;
elseif flags.do_valid
   m = numel(w.g{1}.h);
   a = w.a(1);
   rred = (a^J-1)/(a-1)*(m-a);
   blocksize=w.a(1)^J;
   L=rred+floor((Ls-rred)/blocksize)*blocksize;
else
   L = Ls;
end

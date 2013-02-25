function L=fwtlength(Ls,h,J,varargin);
%FWTLENGTH  FWT length from signal
%   Usage: L=fwtlength(Ls,h,J);
%          L=fwtlength(Ls,h,J,...);
%
%   `fwtlength(Ls,h,J)` returns the length of a Wavelet system that is long
%   enough to expand a signal of length *Ls*. Please see the help on
%   |fwt|_ for an explanation of the parameters *h* and *J*.
%
%   If the returned length is longer than the signal length, the signal
%   will be zero-padded by |fwt|_ to length *L*.
%
%   See also: fwt

% Initialize the wavelet filters structure
h = fwtinit(h,'ana');

definput.import = {'fwt'};
[flags,kv]=ltfatarghelper({},definput,varargin);

if(flags.do_per)
   blocksize=h.a(1)^J;
   L=ceil(Ls/blocksize)*blocksize;
else
   L = Ls;
end

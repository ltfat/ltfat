function g=fir2long(g,Llong);
%FIR2LONG   Extend FIR window to LONG
%   Usage:  g=fir2long(g,Llong);
%
%   FIR2LONG(g,Llong) will extend the FIR window g to a length Llong window
%   by inserting zeros. Note that this is a slightly different behaviour
%   than MIDDLEPAD.
%
%   FIR2LONG can also be used to extend a FIR window to a longer FIR window,
%   for instance in order to satisfy that the window length is divisable by
%   the number of channels.
%
%   See also:  long2fir, middlepad

error(nargchk(2,2,nargin));

Lfir=length(g);


if Lfir>Llong
  error('Llong must be larger than length of window.');
end;

if rem(Lfir,2)==0
  % HPE middlepad works the same way as the FIR extension (e.g. just
  % inserting zeros) for even-length signals.
  g=middlepad(g,Llong,'hp');
else
  % WPE middlepad works the same way as the FIR extension (e.g. just
  % inserting zeros) for odd-length signals.
  g=middlepad(g,Llong);
end;
  

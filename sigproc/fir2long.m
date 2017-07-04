function gout=fir2long(gin,Llong);
%FIR2LONG   Extend FIR window to LONG
%   Usage:  g=fir2long(g,Llong);
%
%   `fir2long(g,Llong)` will extend the FIR window *g* to a length *Llong*
%   window by inserting zeros. Note that this is a slightly different
%   behaviour than |middlepad|.
%
%   `fir2long` can also be used to extend a FIR window to a longer FIR
%   window, for instance in order to satisfy the usual requirement that the
%   window length should be divisible by the number of channels.
%
%   If the input to `fir2long` is a cell, `fir2long` will recurse into
%   the cell array.
%
%   See also:  long2fir, middlepad

complainif_argnonotinrange(nargin,2,2,mfilename);

if iscell(gin)
    gout=cellfun(@(x) fir2long(x,Llong),gin,'UniformOutput',false);    
else
    
    Lfir=length(gin);
    
    if Lfir>Llong
        error('Llong must be larger than length of window.');
    end;
    
    if rem(Lfir,2)==0
        % HPE middlepad works the same way as the FIR extension (e.g. just
        % inserting zeros) for even-length signals.
        gout=middlepad(gin,Llong,'hp');
    else
        % WPE middlepad works the same way as the FIR extension (e.g. just
        % inserting zeros) for odd-length signals.
        gout=middlepad(gin,Llong);
    end;
    
end;


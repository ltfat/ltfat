function g=firfilter(name,M,varargin)
%FIRFILTER  Construct an FIR filter
%   Usage:  g=firfilter(name,M);
%           g=firfilter(name,M,...);
%
%   `firfilter(name,M)` creates an FIR filter of length *M*. This is
%   exactly the same as calling |firwin|. The name must be one of the
%   accepted window types of |firwin|.
%
%   `firfilter(name,M,centre)` constructs a filter with a centre
%   frequency of *centre* measured in normalized frequencies.
%
%   `firfilter` accepts the following optional parameters:
%
%     'fs',fs     If the sampling frequency *fs* is specified then the length
%                 *M* is specified in seconds and the centre frequency
%                 *centre* in Hz.
%
%     'complex'   Make the filter complex valued if the centre frequency
%                 is non-zero. This is the default.
%
%     'real'      Make the filter real-valued if the centre frequency
%                 is non-zero.
%
%     'delay',d   Set the delay of the filter. Default value is zero.
%
%     'causal'    Create a causal filter starting at the first sample. If
%                 specified, this flag overwrites the delay setting.
%
%   It is possible to normalize the impulse response of the filter by
%   passing any of the flags from the |normalize| function. The default
%   normalization is `'area'`, ensuring that the filter has 0dB
%   attenuation at its centre frequency.
%
%   The filter can be used in the |pfilt| routine to filter a signal, or
%   in can be placed in a cell-array for use with |filterbank| or |ufilterbank|.
%
%   See also: blfilter, firwin, pfilt, filterbank

% XXX Implement passing additional parameters to firwin

% Define initial value for flags and key/value pairs.
definput.import={'normalize'};
definput.importdefaults={'area'};
definput.keyvals.delay=0;
definput.keyvals.centre=0;
definput.keyvals.fs=[];
definput.flags.delay={'delay','causal'};
definput.flags.real={'complex','real'};

[flags,kv]=ltfatarghelper({'centre'},definput,varargin);

if ~isempty(kv.fs)
    M=round(M*kv.fs);
    kv.centre=kv.centre/kv.fs*2;
end;

if flags.do_causal
    g.offset=0;
    smallshift=0;
else
    d=floor(kv.delay);
    smallshift=d-floor(d);
    g.offset=d-floor(M/2);
end;

g.h=fftshift(firwin(name,M,'shift',smallshift,flags.norm));
g.centre=kv.centre;
g.realonly=flags.do_real;
g.fs=kv.fs;





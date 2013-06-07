function g=blfilter(name,fsupp,varargin)
%BLFILTER  Construct a band-limited filter
%   Usage:  g=blfilter(name,fsupp,centre);
%           g=blfilter(name,fsupp,centre,...);
%
%   `blfilter(name,fsupp)` construct a band-limited filter. The
%   parameter *name* specifies the shape of the frequency response. The
%   name must be one of the shapes accepted by |firwin|. The support of
%   the frequency response measured in normalized frequencies is
%   specified by *fsupp*.
%
%   `blfilter(name,fsupp,centre)` constructs a filter with a centre
%   frequency of *centre* measured in normalized frequencies.
%
%   `blfilter` accepts the following optional parameters:
%
%     'fs',fs     If the sampling frequency *fs* is specified then the support
%                 *fsupp* and the centre frequency *centre* is specified in Hz.
%
%     'complex'   Make the filter complex valued if the centre frequency
%                 is non-zero.necessary. This is the default.
%
%     'real'      Make the filter real-valued if the centre frequency
%                 is non-zero.
%
%     'delay',d   Set the delay of the filter. Default value is zero.
%
%   It is possible to normalize the transfer function of the filter by
%   passing any of the flags from the |normalize| function. The default
%   normalization is `'peak'`, ensuring that the filter has 0dB
%   attenuation at its centre frequency.
%
%   The filter can be used in the |pfilt| routine to filter a signal, or
%   in can be placed in a cell-array for use with |filterbank| or |ufilterbank|.
%
%   See also: blfilter, firwin, pfilt, filterbank

% Define initial value for flags and key/value pairs.
definput.import={'normalize'};
definput.importdefaults={'peak'};
definput.keyvals.delay=0;
definput.keyvals.centre=0;
definput.keyvals.fs=[];
definput.flags.real={'complex','real'};

[flags,kv]=ltfatarghelper({'centre'},definput,varargin);

if ~isempty(kv.fs)
    fsupp=fsupp/kv.fs*2;
    kv.centre=kv.centre/kv.fs*2;
end;

% Sanitize
kv.centre=modcent(kv.centre,2);

g.H=@(L)    fftshift(firwin(name,round(L/2*fsupp),flags.norm));
g.foff=@(L) round(L/2*kv.centre)-floor(round(L/2*fsupp)/2);
g.realonly=flags.do_real;
g.delay=kv.delay;
g.fs=kv.fs;

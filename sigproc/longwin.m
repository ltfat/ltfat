function h = longwin(name,L,varargin)
%FREQWIN Frequency response window
%   Usage: H = freqwin(name,L,bw);
%
%   `freqwin(name,L,bw)` returns a frequency window *name* of length *L* 
%   with the mainlobe -6dB (half height) bandwidth *bw*. It is intended to
%   represent frequency response of a band-pass filter/window with 
%   bandwidth *bw*. The bandwidth is given in normalised frequencies.
%
%   The function is not periodically wrapped should it be nonzero outside
%   of the *L* samples (as opposed to e.g. |pgauss|).
%
%   The following windows can be used in place of *name*:
%
%     'gauss'        Gaussian window
%
%     'gammatone'    Gammatone window. The default order is 4. Custom order 
%                    can be set by `{'gammatone',order}`.
%
%     'butterworth'  Butterworth window. The default order is 4. Custom 
%                    order can be set by `{'butterworth',order}`.
%
%   `freqwin(name,L,bw,fs)` does the same as above except *bw* is expected
%   to be in Hz given sampling frequency *fs*.
%
%   `freqwin` understands the following key-value pairs and flags at the end of 
%   the list of input parameters:
%
%     'fs',fs      If the sampling frequency *fs* is specified then the *bw*
%                  is expected to be in Hz.
%
%     'shift',s    Shift the window by $s$ samples. The value can be a
%                  fractional number.   
%
%     'wp'         Output is whole point even. This is the default. It
%                  corresponds to a shift of $s=0$.
%
%     'hp'         Output is half point even, as most Matlab filter
%                  routines. This corresponds to a shift of $s=-.5$
%
%   Additionally, the function accepts flags to normalize the output. Please see
%   the help of |normalize|. Default is to use `'peak'` normalization.
%
%
%   See also: firwin, normalize, plotfft

% AUTHORS: Zdenek Prusa

complainif_notenoughargs(nargin,2,upper(mfilename));

if ~isscalar(L)
    error('%s: L must be a scalar',upper(mfilename));
end

longwintypes = arg_longwin(struct);
longwintypes = longwintypes.flags.wintype;

if ~iscell(name), name = {name}; end

if ~ischar(name{1}) || ~any(strcmpi(name{1},longwintypes))
  error(['%s: First input argument must a name of a supported window:',...
         '\n%s'],upper(mfilename), strjoin(longwintypes,'\n'));
end;

winArgs = name(2:end);
winName = lower(name{1});

definput.import={'normalize'};
definput.importdefaults={'null'};
definput.flags.asymmetry = {'asym','sym'};
definput.keyvals.shift = 0;
definput.keyvals.damp = 1;        

switch winName
    case 'reds'
      
    case 'fof'
        definput.keyvals.order = 1;
    case {'damp'}
    otherwise 
        error('%s: SENTINEL. Unknown window.',upper(mfilename));
end

[flags,kv]=ltfatarghelper({},definput,varargin);

if kv.damp < 0
    error('%s: damp must be nonnegative.',upper(mfilename));
end

switch winName
    case 'reds'
      
    case 'fof'
end

h = zeros(L,1);
switch winName
    case 'damp'
        l = 0:floor(L/2);
        h(1:numel(l)) = exp(-kv.damp*l);
    case 'reds'
        

    case 'fof'
        delay = (1/kv.damp)*acos();
        l = 0:L-1;
        h = 0.5*(1-cos(kv.damp*l));
end

h=normalize(h,flags.norm);


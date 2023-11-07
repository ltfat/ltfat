function [g, info] = pwin(name,varargin)
%PWIN periodically sampled window
%
%   Input parameters:
%     name   : window name
%     tshift : time shift
%     fshift : frequency shift
%
%   Output parameters:
%     g    : window
%     info : struct with additional parameters
%
%   `pwin` generates a periodically sampled cosine-type windows that can be
%   fractionally modulated and shifted. Per default, the window is designed 
%   such that one period of the underlying cosine fits into its nominal
%   length `L`. `fshift` decreases that period additively. The `tshift`
%   indicates the factor times pi by which the underlying cosine is shifted.
%   Currently supported windows: 'hanning', 'hamming'
%
%   Optional parameters:
%     'tshift',tshift   shift in time ([0,1]...half a period)
%     'fshift',fshift   shift in frequency (max depends on L)
%     'fs',fs           sampling frequency, controls shift resolution
%
%   See also:  freqwin, pgauss, firwin, firkaiser


definput.import={'setnorm'};
definput.importdefaults={'null'};
definput.flags.centering={'wp','hp','shift'};
definput.keyvals.fs = 1000;
definput.keyvals.fc = 1;
definput.keyvals.phi = 0;
definput.keyvals.L = 1;
definput.keyvals.tshift = 0;
definput.keyvals.fshift = 0;


[flags,kv]=ltfatarghelper({},definput,varargin);


if ischar(name)
    kv.phi = -2*pi*kv.fc*kv.tshift;
    kv.fc = kv.fc + kv.fshift;

    %make the cosine...
    dt = 1/kv.fs; %the cosine should always have length T*kv.fs
    t = (0:dt:kv.L-dt)';
    x = cos(2*pi*kv.fc*t + kv.phi) + 1; % '+1' is just the offset, quasi-arbitrary

    %...make the rectangular window...
    T = 1/kv.fc;
    rect = zeros(size(x));
    rect(1:floor(T*kv.fs)) = ones(floor(T*kv.fs),1) * 0.5;%this '0.5' is the counterpart to the offset above
    rect = circshift(rect, round(kv.tshift*T*kv.fs));
    %...and the periodically sampled signal
    persamp = x.*rect;
    persamp = middlepad(persamp(persamp~=0), numel(persamp));
    persamp = circshift(persamp, round(kv.tshift*T*kv.fs));

    switch name  
     %the windows are then more or less their superpositions
     case {'hanning'}
         g=(0.5+0.5*persamp);
     case 'hamming'
         g=0.54+0.46*persamp;
     otherwise
      error('Unknown window: %s.',name);
    end
else
    g = [];
    disp('not yet implemented.');
end
    g = setnorm(g, flags.norm);
      
info = [];
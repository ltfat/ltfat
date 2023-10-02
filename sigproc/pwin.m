function [g, info] = pwin(name,tshift,fshift,varargin)
%PWIN periodically sampled window
%
%   Input parameters:
%     name : window name
%     tshift : time shift
%     fshift : frequency shift
%
%   Output parameters:
%     g : window
%


definput.import={'setnorm'};
definput.keyvals.fs = 10000;
definput.keyvals.fc = 10;
definput.keyvals.phi = 0;
definput.keyvals.L = 1;


[flags,kv]=ltfatarghelper({},definput,varargin);

kv.phi = 2*pi*kv.fc*tshift;
kv.fc = kv.fc + fshift;

%make the cosine...
dt = 1/kv.fs;
t = (0:dt:kv.L-dt)';
x = cos(2*pi*kv.fc*t + kv.phi) + 1; % '+1' is just the offset, quasi-arbitrary

%...make the rectangular window...
T = 1/kv.fc;
rect = zeros(size(x));
rect(1:floor(T*kv.fs)) = ones(floor(T*kv.fs),1) * 0.5;%this '0.5' is the counterpart to the offset above
rect = circshift(rect, floor(kv.phi/(2*pi*kv.fc)*kv.fs));

%...and the periodically sampled signal
persamp = x.*rect;

switch name  
 %the windows are then more or less their superpositions
 case {'hanning'}
     g=(0.5+0.5*persamp);
 otherwise
  error('Unknown window: %s.',name);
end;

g=setnorm(g,flags.norm);

info = [];
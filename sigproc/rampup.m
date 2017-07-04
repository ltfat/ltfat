function outsig=rampup(L,wintype)
%RAMPUP  Rising ramp function
%   Usage: outsig=rampup(L);
%
%   `rampup(L)` will return a rising ramp function of length *L*. The
%   ramp is a sinusoide starting from zero and ending at one. The ramp
%   is centered such that the first element is always 0 and the last
%   element is not quite 1, such that the ramp fits with following ones.
%
%   `rampup(L,wintype)` will use another window for ramping. This may be any
%   of the window types from |firwin|. Please see the help on |firwin| for
%   more information. The default is to use a piece of the Hann window.
%
%   See also: rampdown, rampsignal, firwin

complainif_argnonotinrange(nargin,1,2,mfilename);

if nargin==1
  wintype='hann';
end;
  
win=firwin(wintype,2*L,'inf');
outsig=win(L+1:2*L);


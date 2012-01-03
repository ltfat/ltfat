function outsig=rampup(L,wintype)
%RAMPUP  Rising ramp function
%   Usage: outsig=rampup(L);
%
%   RAMPUP(L) will return a rising ramp function of length L. The
%   ramp is a sinusoide starting from zero and ending at one. The ramp
%   is centered such that the first element is always 0 and the last
%   element is not quite 1, such that the ramp fits with following ones.
%
%   RAMPUP(L,wintype) will use another window for ramping. This may be
%   any of the window types from FIRWIN. Please see the help on FIRWIN
%   for more information. The default is to use a piece of the Hann window.
%
%   See also: rampdown, rampsignal, firwin

error(nargchk(1,2,nargin))

if nargin==1
  wintype='hann';
end;
  
win=firwin(wintype,2*L,'inf');
outsig=win(L+1:2*L);
  
%OLDFORMAT

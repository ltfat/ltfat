function outsig=rampdown(L,wintype)
%RAMPDOWN  Falling ramp function
%   Usage: outsig=rampdown(siglen);
%
%   `rampdown(siglen)` will return a falling ramp function of length
%   *siglen*. The ramp is a sinusoid starting from one and ending at
%   zero. The ramp is centered such that the first element is always one and
%   the last element is not quite zero, such that the ramp fits with
%   following zeros.
%
%   `rampdown(L,wintype)` will use another window for ramping. This may be
%   any of the window types from |firwin|. Please see the help on |firwin|
%   for more information. The default is to use a piece of the Hann window.
%  
%   See also: rampup, rampsignal, firwin
   
complainif_argnonotinrange(nargin,1,2,mfilename);

if nargin==1
  wintype='hann';
end;
  
win=firwin(wintype,2*L,'inf');
outsig=win(1:L);


function freq = audtofreq(aud,varargin)
%AUDTOFREQ  Converts auditory units to frequency (Hz)
%   Usage: freq = audtofreq(aud);
%  
%   `audtofreq(aud,scale)` converts values on the selected auditory scale to
%   frequencies measured in Hz.
%
%   See the help on |freqtoaud| to get a list of the supported values of
%   the *scale* parameter. If no scale is given, the erb-scale will be
%   selected by default.
%
%   See also: freqtoaud, audspace, audfiltbw

%   AUTHOR: Peter L. Søndergaard

%% ------ Checking of input parameters ---------

if nargin<1
  error('%s: Too few input parameters.',upper(mfilename));
end;


if ~isnumeric(aud) 
  error('%s: aud must be a number.',upper(mfilename));
end;

definput.import={'freqtoaud'};
[flags,kv]=ltfatarghelper({},definput,varargin);

%% ------ Computation --------------------------
  
if flags.do_mel
  freq = 700*sign(aud).*(exp(abs(aud)*log(17/7)/1000)-1);
end;

if flags.do_mel1000
  freq = 1000*sign(aud).*(exp(abs(aud)*log(2)/1000)-1);
end;

if flags.do_erb
  freq = (1/0.00437)*sign(aud).*(exp(abs(aud)/9.2645)-1);
end;

if flags.do_bark
  % This one was found through http://www.ling.su.se/STAFF/hartmut/bark.htm
  freq = sign(aud).*1960./(26.81./(abs(aud)+0.53)-1);
end;

if flags.do_erb83
  %freq = 14363./(1-exp((aud-43.0)/11.7))-14675;
  freq = sign(aud).*14675.*(1-exp(0.0895255*abs(aud)))./(exp(0.0895255*abs(aud)) - 47.0353);
end;

if flags.do_freq
  freq=aud;
end;

if flags.do_log10 || flags.do_semitone
   freq = 10^(aud);
end

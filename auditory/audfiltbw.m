function bw = audfiltbw(fc,varargin)
%AUDFILTBW  Bandwidth of auditory filter
%   Usage: bw = audfiltbw(fc)
%
%   `audfiltbw(fc)` returns the critical bandwidth of the auditory filter 
%   at center frequency *fc* defined in equivalent rectangular bandwidth.
%   The function uses the relation:
%
%   .. bw = 24.7 + fc/9.265
%
%   .. math::  bw = 24.7 + \frac{fc}{9.265}
%     
%   as estimated in Glasberg and Moore (1990). This function is also used
%   when the original ERB scale ('erb83') is chosen.
% 
%   `audfiltbw(fc,'bark')` returns the critical bandwidth at *fc* according
%   to the Bark scale using the relation:
% 
%   .. bw = 25 + 75 ( 1+1.4*10^{-6} fc^2 )^0.69
%
%   .. math::  bw = 25 + 75 ( 1+1.4\times10^{-6} fc^2 )^{0.69}
% 
%   as estimated by Zwicker and Terhardt (1980).
%
%   For the scales 'mel', 'mel1000', 'log10' and 'semitone', no critical
%   bandwidth function is usually given. Following the example of the
%   equivalent rectangular bandwidth (associated with the ERB scale), we
%   use the derivative of the inverse of the scale function F_{Scale},
%   evaluated at F_{Scale}(fc), i.e. 
% 
%   .. bw = (F_{scale}^{-1})'(F_{scale}(fc))
%
%   .. math::  bw = (F_{scale}^{-1})'(F_{scale}(fc))
%
%   See also: freqtoerb, erbspace, freqtoaud, audtofreq
%
%   References: glasberg1990daf zwickerterhardt80
  
%   AUTHOR : Peter L. Søndergaard
  
% ------ Checking of input parameters ---------
  
% narginchk(1,1);
if nargin<1
  error('%s: Too few input parameters.',upper(mfilename));
end;

if ~isnumeric(fc) || any(fc(:)<0)
  error('AUDFILTBW: fc must be non-negative.');
end;

definput.flags.audscale={'erb','erb83','bark','mel','mel1000','log10','semitone'};
[flags,kv]=ltfatarghelper({},definput,varargin);

% ------ Computation --------------------------

% FIXME: What is the upper frequency for which the estimation is valid?
if flags.do_erb || flags.do_erb83
    bw = 24.7 + fc/9.265;
end

if flags.do_bark
    bw = 25 + 75*(1 + 1.4E-6*fc.^2).^0.69;
end

if flags.do_mel    
    bw = log(17/7)*(700+fc)/1000;
end

if flags.do_mel1000  
    bw = log(2)*(1+fc/1000);
end;

if flags.do_log10 || flags.do_semitone
    bw = fc;
end

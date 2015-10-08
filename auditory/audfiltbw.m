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
%   as estimated in Glasberg and Moore (1990).
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
%   See also: freqtoerb, erbspace
%
%   References: glasberg1990daf zwickerterhardt80
  
%   AUTHOR : Peter L. Søndergaard
  
% ------ Checking of input parameters ---------
  
% error(nargchk(1,1,nargin));
if nargin<1
  error('%s: Too few input parameters.',upper(mfilename));
end;

if ~isnumeric(fc) || any(fc(:)<0)
  error('AUDFILTBW: fc must be non-negative.');
end;

definput.flags.audscale={'erb','bark'};
[flags,kv]=ltfatarghelper({},definput,varargin);

% ------ Computation --------------------------

% FIXME: What is the upper frequency for which the estimation is valid?
if flags.do_erb
    bw = 24.7 + fc/9.265;
end

if flags.do_bark
    bw = 25 + 75*(1 + 1.4E-6*fc.^2).^0.69;
end

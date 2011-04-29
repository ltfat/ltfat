function [y,n] = audspacebw(flow,fhigh,varargin)
%AUDSPACEBW  Auditory scale points specified by bandwidth
%   Usage: y=audspacebw(flow,fhigh,bw,hitme);
%          y=audspacebw(flow,fhigh,bw);
%          y=audspacebw(flow,fhigh);
%          [y,n]=audspacebw(...);
%
%   AUDSPACEBW(flow,fhigh,bw,scale) computes a vector containing
%   values equistantly scaled between flow and fhigh on the selected
%   auditory scale.  All frequencies are specified in Hz.The distance
%   between two consecutive values are given by a value of bw the
%   scale, and the points will be centered on the scale between flow
%   and fhigh.
%
%   See the help on FREQTOAUD to get a list of the supported values of the
%   scale parameter.
%  
%   AUDSPACEBW(flow,fhigh,bw,hitme,scale) will do as above, but one of
%   the points is quaranteed to be the frequency hitme.
%
%   [y,n]=AUDSPACEBW( ... ) additionally returns the number of points n in
%   the output vector y.
%
%   See also: freqtoaud, audspace, audfiltbw
  
%   AUTHOR : Peter L. Soendergaard
  
% ------ Checking of input parameters ---------
  
if nargin<2
  error('%s: Too few input parameters.',upper(mfilename));
end;

if ~isnumeric(flow) || ~isscalar(flow) || flow<0
  error('%s: flow must be a non-negative scalar.',upper(mfilename));
end;

if ~isnumeric(fhigh) || ~isscalar(fhigh) || fhigh<0
  error('%s: fhigh must be a non-negative scalar.',upper(mfilename));
end;

if flow>fhigh
  error('%s: flow must be less than or equal to fhigh.',upper(mfilename));
end;

definput.import={'freqtoaud'};
definput.keyvals.hitme=[];
definput.keyvals.bw=1;

[flags,kv,bw]=ltfatarghelper({'bw','hitme'},definput,varargin);

if ~isnumeric(bw) || ~isscalar(bw) || bw<=0 
  error('%s: bw must be a positive scalar.',upper(mfilename));
end;

  
%% ------ Computation --------------------------

if isempty(kv.hitme)
  % Convert the frequency limits to auds.
  audlimits = freqtoaud([flow,fhigh],flags.audscale);
  audrange  = audlimits(2)-audlimits(1);

  % Calculate number of points, excluding final point
  n         = floor(audrange/bw);

  % The remainder is calculated in order to center the points
  % correctly between flow and fhigh.
  remainder = audrange-n*bw;

  audpoints = audlimits(1)+(0:n)*bw+remainder/2;
  
  % Add the final point
  n=n+1;  
  
else
    
  % Convert the frequency limits to auds.
  audlimits    = freqtoaud([flow,fhigh,kv.hitme],flags.audscale);
  audrangelow  = audlimits(3)-audlimits(1);
  audrangehigh = audlimits(2)-audlimits(3);

  % Calculate number of points, exluding final point.
  nlow = floor(audrangelow/bw);
  nhigh = floor(audrangehigh/bw);
  
  audpoints=(-nlow:nhigh)*bw+audlimits(3);
  n=nlow+nhigh+1;
end;

y = audtofreq(audpoints,flags.audscale);

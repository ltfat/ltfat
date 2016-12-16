function [y,n] = audspacebw(fmin,fmax,varargin)
%AUDSPACEBW  Auditory scale points specified by bandwidth
%   Usage: y=audspacebw(fmin,fmax,bw,hitme);
%          y=audspacebw(fmin,fmax,bw);
%          y=audspacebw(fmin,fmax);
%          [y,n]=audspacebw(...);
%
%   `audspacebw(fmin,fmax,bw,scale)` computes a vector containing values
%   equistantly scaled between frequencies *fmin* and *fmax* on the
%   selected auditory scale.  All frequencies are specified in Hz.The
%   distance between two consecutive values is *bw* on the selected scale,
%   and the points will be centered on the scale between *fmin* and *fmax*.
%
%   See the help on |freqtoaud| to get a list of the supported values of the
%   *scale* parameter.
%  
%   `audspacebw(fmin,fmax,bw,hitme,scale)` will do as above, but one of
%   the points is quaranteed to be the frequency *hitme*.
%
%   `[y,n]=audspacebw(...)` additionally returns the number of points *n* in
%   the output vector *y*.
%
%   See also: freqtoaud, audspace, audfiltbw
  
%   AUTHOR : Peter L. Søndergaard
  
% ------ Checking of input parameters ---------
  
if nargin<2
  error('%s: Too few input parameters.',upper(mfilename));
end;

if ~isnumeric(fmin) || ~isscalar(fmin) || fmin<0
  error('%s: fmin must be a non-negative scalar.',upper(mfilename));
end;

if ~isnumeric(fmax) || ~isscalar(fmax) || fmax<0
  error('%s: fmax must be a non-negative scalar.',upper(mfilename));
end;

if fmin>fmax
  error('%s: fmin must be less than or equal to fmax.',upper(mfilename));
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
  audlimits = freqtoaud([fmin,fmax],flags.audscale);
  audrange  = audlimits(2)-audlimits(1);

  % Calculate number of points, excluding final point
  n         = floor(audrange/bw);

  % The remainder is calculated in order to center the points
  % correctly between fmin and fmax.
  remainder = audrange-n*bw;

  audpoints = audlimits(1)+(0:n)*bw+remainder/2;
  
  % Add the final point
  n=n+1;  
  
else
    
  % Convert the frequency limits to auds.
  audlimits    = freqtoaud([fmin,fmax,kv.hitme],flags.audscale);
  audrangelow  = audlimits(3)-audlimits(1);
  audrangehigh = audlimits(2)-audlimits(3);

  % Calculate number of points, exluding final point.
  nlow = floor(audrangelow/bw);
  nhigh = floor(audrangehigh/bw);
  
  audpoints=(-nlow:nhigh)*bw+audlimits(3);
  n=nlow+nhigh+1;
end;

y = audtofreq(audpoints,flags.audscale);

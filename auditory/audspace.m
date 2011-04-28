function [y,bw] = audspace(flow,fhigh,n,varargin)
%AUDSPACE  Equidistantly spaced points on auditory scale
%   Usage: y=audspace(scale,flow,fhigh,n);
%
%   AUDSPACE(flow,fhigh,n,scale) computes a vector of length n containing values
%   equistantly scaled on the selected auditory scale between flow and fhigh. All
%   frequencies are specified in Hz.
%
%   See the help on FREQTOAUD to get a list of the supported values of the
%   scale parameter.
%  
%   [y,bw]=AUDSPACE(...) does the same but outputs the bandwidth between
%   each sample measured on the selected scale.
%  
%   See also: freqtoaud, audspacebw, audfiltbw
%
%R  glasberg1990daf
  
%   AUTHOR : Peter L. Soendergaard
  
%% ------ Checking of input parameters ---------

if nargin<3
  error('%s: Too few input parameters.',upper(mfilename));
end;
  
% Default parameters.

if ~isnumeric(flow) || ~isscalar(flow) || flow<0
  error('%s: flow must be a non-negative scalar.',upper(mfilename));
end;

if ~isnumeric(fhigh) || ~isscalar(fhigh) || fhigh<0
  error('%s: fhigh must be a non-negative scalar.',upper(mfilename));
end;

if ~isnumeric(n) || ~isscalar(n) || n<=0 || fix(n)~=n
  error('%s: n must be a positive, integer scalar.',upper(mfilename));
end;

if flow>fhigh
  error('%s: flow must be less than or equal to fhigh.',upper(mfilename));
end;

definput.import={'freqtoaud'};
[flags,kv]=ltfatarghelper({},definput,varargin);


%% ------ Computation --------------------------

audlimits = freqtoaud([flow,fhigh],flags.audscale);

y = audtofreq(linspace(audlimits(1),audlimits(2),n),flags.audscale);

bw=(audlimits(2)-audlimits(1))/(n-1);

% Set the endpoints to be exactly what the user specified, instead of the
% calculated values
y(1)=flow;
y(end)=fhigh;
function [] = plotnsdgtreal(coef,a,varargin)
%PLOTNSDGTREAL Plot spectrogram from NSDGTREAL coefficients
%   Usage:  plotnsdgtreal(c,a,dynrange,sr);
%
%   Input parameters:
%         coef     : Cell array of coefficients.
%         a        : Vector of time positions of windows.
%         dynrange : Colorscale dynamic range in dB (default 60 dB).
%         sr       : signal sample rate in Hz (default 1 Hz).
%
%   PLOTNSDGTREAL(coef,a) plots coefficients computed NSDGTREAL. For
%   more details on the format of the variables coef and _a, please read the
%   NSDGT function help.
%
%   PLOTNSDGTREAL(coef,a,fs) will do the same assuming a sampling rate of
%   fs Hz of the original signal.
%
%   PLOTNSDGTREAL(coef,a,fs,dynrange) will additionally limit the dynamic
%   range.
%
%   PLOTNSDGTREAL supports all the optional parameters of TFPLOT. Please
%   see the help of TFPLOT for an exhaustive list. In addition, the
%   following parameters may be specified:
%
%-  'xres',xres - Approximate number of pixels along x-axis /time.
%               Default value is 800
%
%-  'yres',yres - Approximate number of pixels along y-axis / frequency
%               Default value is 600
%
%   See also: tfplot, nsdgt, nsdgtreal

%   AUTHOR : Florent Jaillet & Peter L. Soendergaard
%   TESTING: OK 
%   REFERENCE: NA

if nargin<2
  error('%s: Too few input parameters.',upper(mfilename));
end;

definput.import={'ltfattranslate','tfplot'};

definput.keyvals.xres=800;
definput.keyvals.yres=600;

[flags,kv,fs]=ltfatarghelper({'fs','dynrange'},definput,varargin);

timepos=cumsum(a)-a(1);

N=length(a);
cwork=zeros(kv.yres,N);

%% -------- Interpolate in frequency ---------------------

for ii=1:N
  column=coef{ii};
  M=length(column);
  cwork(:,ii)=interp1(linspace(0,1,M),column,linspace(0,1,kv.yres),'nearest');
end;

%% --------  Interpolate in time -------------------------

% Time step in next equidistant spacing on the x-axis (in samples)
aplot=timepos(end)/kv.xres;

% Time positions where we want our pixels plotted (in samples)
xr=(0:kv.xres-1)*aplot;

coef=zeros(kv.yres,kv.xres);
for ii=1:kv.yres
  data=interp1(timepos,cwork(ii,:).',xr,'nearest').';
  coef(ii,:)=data;
end;

yr=[0,1];

tfplot(coef,aplot,yr,'argimport',flags,kv);

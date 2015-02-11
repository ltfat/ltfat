function coef = plotnsdgtreal(coef,a,varargin)
%PLOTNSDGTREAL Plot NSDGTREAL coefficients
%   Usage:  plotnsdgtreal(c,a,fs,dynrange);
%
%   Input parameters:
%         coef     : Cell array of coefficients.
%         a        : Vector of time positions of windows.
%         fs       : signal sample rate in Hz (optional).
%         dynrange : Colorscale dynamic range in dB (optional).
%
%   `plotnsdgtreal(coef,a)` plots coefficients computed using |nsdgtreal| or
%   |unsdgtreal|. For more details on the format of the variables *coef* and *a*,
%   please read the function help for these functions.
%
%   `plotnsdgtreal(coef,a,fs)` does the same assuming a sampling rate of
%   *fs* Hz of the original signal.
%
%   `plotnsdgtreal(coef,a,fs,dynrange)` additionally limits the dynamic range.
%
%   `C=plotnsdgtreal(...)` returns the processed image data used in the
%   plotting. Inputting this data directly to `imagesc` or similar
%   functions will create the plot. This is useful for custom
%   post-processing of the image data.
%
%   `plotnsdgtreal` supports all the optional parameters of |tfplot|. Please
%   see the help of |tfplot| for an exhaustive list. In addition, the
%   following parameters may be specified:
%
%     'xres',xres  Approximate number of pixels along x-axis /time.
%                  Default value is 800
%
%     'yres',yres  Approximate number of pixels along y-axis / frequency
%                  Default value is 600
%
%   See also: tfplot, nsdgt, nsdgtreal

%   AUTHOR : Florent Jaillet & Peter L. Søndergaard
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

coef=tfplot(coef,aplot,yr,'argimport',flags,kv);

if nargout<1
    clear coef;
end

function [] = plotnsdgt(c,a,fs,varargin)
%PLOTNSDGT Plot spectrogram from nonstationary Gabor coefficients
%   Usage:  plotnsdgt(c,a,dynrange,sr);
%
%   Input parameters:
%         c        : Cell array of coefficients.
%         a        : Vector of time positions of windows.
%         dynrange : Colorscale dynamic range in dB (default 60 dB).
%         sr       : signal sample rate in Hz (default 1 Hz).
%
%   PLOTNSDGT the spectrogram from coefficients computed with the functions
%   NSDGT or NSDGTREAL. For more details on the format of the variables c
%   and _a format, please read the NSDGT function help.
%
%   Additional arguments can be supplied like this:
%   NSSGRAM(f,'nf','tfr',2,'log'). The arguments must be character strings
%   possibly followed by an argument:
%
%-  'real'    - Assume coefficients from NSDGTREAL. This is the default.
%
%-  'complex' - Assume coefficients from NSDGT.
%
%-  'image'   - Use 'imagesc' to display the spectrogram. This is the
%               default.
%
%-  'clim',[clow,chigh] - Use a colormap ranging from clow to chigh. These
%               values are passed to IMAGESC. See the help on IMAGESC.
%
%-  'dynrange',r - Use a colormap in the interval [chigh-r,chigh], where
%               chigh is the highest value in the plot.
%
%-  'xres',xres - Approximate number of pixels along x-axis / time.
%
%-  'yres',yres - Approximate number of pixels along y-axis / frequency
%
%-  'contour' - Do a contour plot to display the spectrogram.
%          
%-  'surf'    - Do a surf plot to display the spectrogram.
%
%-  'mesh'    - Do a mesh plot to display the spectrogram.
%
%-  'colorbar' - Display the colorbar. This is the default.
%
%-  'nocolorbar' - Do not display the colorbar.
%
%   See also: nsdgt, nsdgtreal, nssgram

%   AUTHOR : Florent Jaillet & Peter L. Soendergaard
%   TESTING: OK 
%   REFERENCE: NA

timepos=cumsum(a)-a(1);

% Define initial value for flags and key/value pairs.
defnopos.flags.plottype={'image','contour','mesh','pcolor'};
defnopos.flags.transformtype={'real','complex'};

defnopos.flags.clim={'noclim','clim'};
defnopos.flags.colorbar={'colorbar','nocolorbar'};

defnopos.keyvals.clim=[0,1];
defnopos.keyvals.dynrange=[];
defnopos.keyvals.xres=800;
defnopos.keyvals.yres=600;

[flags,kv]=ltfatarghelper({'dynrange'},defnopos,varargin,mfilename);

cwork=zeros(kv.yres,length(a));

%% -------- Interpolate in frequency ---------------------

for ii=1:length(a)
  column=20*log10(abs(c{ii}+realmin));
  M=length(column);
  cwork(:,ii)=interp1(linspace(0,1,M),column,linspace(0,1,kv.yres),'nearest');
end;

%% --------  Interpolate in time -------------------------
% this is non-equidistant, so we use a cubic spline

% Time positions (in Hz) of our samples.
timepos = (cumsum(a)-a(1))/fs;

% Time positions where we want our pixels plotted.
xr=((0:kv.xres-1)/kv.xres*timepos(end)).';

coef=zeros(kv.yres,kv.xres);
for ii=1:kv.yres
  data=interp1(timepos,cwork(ii,:).',xr,'nearest').';
  coef(ii,:)=data;
end;

% 'dynrange' parameter is handled by thresholding the coefficients.
if ~isempty(kv.dynrange)
  maxclim=max(coef(:));
  coef(coef<maxclim-kv.dynrange)=maxclim-kv.dynrange;
end;

xr=linspace(0,timepos(end),kv.xres);
yr=linspace(0,fs/2,kv.yres);
  
switch(flags.plottype)
  case 'image'
    if flags.do_clim
      imagesc(xr,yr,coef,kv.clim);
    else
     imagesc(xr,yr,coef);
    end;
  case 'contour'
    contour(xr,yr,coef);
  case 'surf'
    surf(xr,yr,coef);
  case 'pcolor'
    pcolor(xr,yr,coef);
end;

if flags.do_colorbar
  colorbar;
end;

axis('xy');
xlabel('Time (s)')
ylabel('Frequency (Hz)')



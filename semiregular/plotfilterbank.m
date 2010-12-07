function [] = plotfilterbank(coef,a,varargin)
%PLOTNSDGT Plot spectrogram from filterbank coefficients
%   Usage:  plotnsdgt(c,a,dynrange,sr);
%
%   Input parameters:
%         c        : Cell array of coefficients.
%         a        : Vector of time positions of windows.
%         dynrange : Colorscale dynamic range in dB (default 60 dB).
%         sr       : signal sample rate in Hz (default 1 Hz).
%
%   PLOTNSDGT(c,a) plots the spectrogram from coefficients computed with the
%   function FILTERBANK. For more details on the format of the variables c
%   and _a, please read the NSDGT function help.
%
%   The function takes the following arguments at the end of the command line:
%
%     'fs'         - Assume a sampling rate of fs Hz.
%%
%-    'image'      - Use 'imagesc' to display the spectrogram. This is the
%                    default.
%
%-    'clim',[clow,chigh] - Use a colormap ranging from clow to chigh. These
%                    values are passed to IMAGESC. See the help on IMAGESC.
%
%-    'dynrange',r - Use a colormap in the interval [chigh-r,chigh], where
%                    chigh is the highest value in the plot.
%
%-    'xres',xres  - Approximate number of pixels along x-axis / time.
%
%-    'yres',yres  - Approximate number of pixels along y-axis / frequency
%
%-    'contour'    - Do a contour plot to display the spectrogram.
%          
%-    'surf'       - Do a surf plot to display the spectrogram.
%
%-    'mesh'       - Do a mesh plot to display the spectrogram.
%
%-    'colorbar'   - Display the colorbar. This is the default.
%
%-    'nocolorbar' - Do not display the colorbar.
%
%   See also: filterbank, nssgram

%   AUTHOR : Peter L. Soendergaard
%   TESTING: OK 
%   REFERENCE: NA

% Define initial value for flags and key/value pairs.
definput.flags.plottype={'image','contour','mesh','pcolor'};

definput.flags.clim={'noclim','clim'};
definput.flags.colorbar={'colorbar','nocolorbar'};

definput.keyvals.clim=[0,1];
definput.keyvals.dynrange=[];
definput.keyvals.xres=800;
definput.keyvals.yres=600;
definput.keyvals.fs=[];

[flags,kv,fs]=ltfatarghelper({'fs','dynrange'},definput,varargin);

N=size(coef,1);
M=size(coef,2);

coef=20*log10(abs(coef));

if isempty(fs)
  fs=1;
end;

% Time positions (in Hz) of our samples.
xr = (0:a:N*a-1)/fs;

% Freq. pos is just number of the channel.
yr=1:M;

% 'dynrange' parameter is handled by thresholding the coefficients.
if ~isempty(kv.dynrange)
  maxclim=max(coef(:));
  coef(coef<maxclim-kv.dynrange)=maxclim-kv.dynrange;
end;

  
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



function nssgram(c,a,M,fs,varargin)
%NSSGRAM  Spectrogram using non-stationary Gabor frame.
%   Usage: nssgram(f,op1,op2, ... );
%          nssgram(f,fs,op1,op2, ... );
%          C=nssgram(f, ... );
%
%   NSSGRAM() plots a spectrogram of f using a DGT.
%
%   NSSGRAM(f,fs) does the same for a signal with sampling rate fs (sampled
%   with fs samples per second);
%
%   C=NSSGRAM(f, ... ) returns the image to be displayed as a matrix. Use
%   this in conjunction with IMWRITE etc. Do NOT use this as a method to
%   create a Gabor transform, use DGT or DGTREAL for that instead.
%
%   Additional arguments can be supplied like this:
%   NSSGRAM(f,'nf','tfr',2,'log'). The arguments must be character strings
%   possibly followed by an argument:
%
%-  'db'      - Apply 20*log10 to the coefficients. This makes it possible to
%               see very weak phenomena, but it might show too much noise. A
%               logarithmic scale is more adapted to perception of sound.
%               This is the default.
%
%-  'lin'     - Show the energy of the coefficients on a linear scale.
%
%-  'tfr',v   - Set the ratio of frequency resolution to time resolution.
%               A value v=1 is the default. Setting v>1 will give better
%               frequency resolution at the expense of a worse time
%               resolution. A value of 0<v<1 will do the opposite.
%
%-  'wlen',s  - Window length. Specifies the length of the window
%               measured in samples. See help of PGAUSS on the exact
%               details of the window length.
%
%-  'nf'      - Display negative frequencies, with the zero-frequency
%               centered in the middle. For real signals, this will just
%               mirror the upper half plane. This is standard for complex
%               signals.
%
%-  'tc'      - Time centering. Move the beginning of the signal to the
%               middle of the plot. This is useful for visualizing the
%               window functions of the toolbox.
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
%-  'thr',r   - Keep only the largest fraction r of the coefficients, and
%               set the rest to zero.
%
%-  'fmax',y  - Display y as the highest frequency.
%
%-  'xres',xres - Approximate number of pixels along x-axis / time.
%
%-  'yres',yres - Approximate number of pixels along y-axis / frequency
%                 If only one of 'xres' and
%                 'yres' is specified, the default aspect ratio will be
%                 used.
%
%-  'contour' - Do a contour plot to display the spectrogram.
%          
%-  'surf'    - Do a surf plot to display the spectrogram.
%
%-  'mesh'    - Do a mesh plot to display the spectrogram.
%
  %   See also:  dgt, dgtreal

%   AUTHOR : Peter Soendergaard.
%   TESTING: NA
%   REFERENCE: NA
  
if nargin<1
  error('Too few input arguments.');
end;

% Define initial value for flags and key/value pairs.
defnopos.flags.thr={'nothr','thr'};
defnopos.flags.plottype={'image','contour','mesh','pcolor'};

defnopos.flags.clim={'noclim','clim'};
defnopos.flags.log={'db','lin'};
defnopos.flags.dynrange={'nodynrange','dynrange'};
defnopos.flags.colorbar={'colorbar','nocolorbar'};

defnopos.keyvals.fs=[];
defnopos.keyvals.tfr=1;
defnopos.keyvals.wlen=0;
defnopos.keyvals.thr=0;
defnopos.keyvals.clim=[0,1];
defnopos.keyvals.climsym=1;
defnopos.keyvals.fmax=0;
defnopos.keyvals.dynrange=100;
defnopos.keyvals.xres=800;
defnopos.keyvals.yres=600;

[flags,keyvals]=ltfatarghelper({},defnopos,varargin,mfilename);

cwork=zeros(keyvals.yres,length(a));

%% -------- Interpolate in frequency ---------------------
%  this is equidistant, so we use dctresample 
for ii=1:length(a)
  cwork(:,ii)=dctresample(20*log10(abs(c{ii}+realmin)),keyvals.yres);
end;

%% --------  Interpolate in time -------------------------
% this is non-equidistant, so we use a cubic spline

% Time positions (in Hz) of our samples.
timepos = (cumsum(a)-a(1))/fs;

% Time positions where we want our pixels plotted.
xr=((0:keyvals.xres-1)/keyvals.xres*timepos(end)).';
yr=linspace(0,fs/2,keyvals.yres);

coef=zeros(keyvals.yres,keyvals.xres);
for ii=1:keyvals.yres
  %PP=spline(timepos,cwork(ii,:).');  
  PP=interp1(timepos,cwork(ii,:).','linear','pp');  
  data=ppval(PP,xr).';
  coef(ii,:)=data;
end;
 
%coef=coef(:,40:end-40);

% 'dynrange' parameter is handled by thresholding the coefficients.
if isfield(flags,'dynrange') && flags.do_dynrange
  maxclim=max(coef(:));
  coef(coef<maxclim-keyvals.dynrange)=maxclim-keyvals.dynrange;
end;

if 0
switch(flags.plottype)
  case 'image'
    if flags.do_clim
      imagesc(xr,yr,coef,keyvals.clim);
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
end;

imagesc(coef);

if flags.do_colorbar
  colorbar;
end;

axis('xy');
xlabel('Time (s)')
ylabel('Frequency (Hz)')

function tfplot(coef,a,M,L,resamp,keyvals,flags)
%SGRAM  Spectrogram.
%   Usage: sgram(f,op1,op2, ... );
%          sgram(f,fs,op1,op2, ... );
%          C=sgram(f, ... );
%
%   SGRAM(f) plots a spectrogram of f using a DGT.
%
%   SGRAM(f,fs) does the same for a signal with sampling rate fs (sampled
%   with fs samples per second);
%
%   C=SGRAM(f, ... ) returns the image to be displayed as a matrix. Use
%   this in conjunction with IMWRITE etc. Do NOT use this as a method to
%   create a Gabor transform, use DGT or DGTREAL for that instead.
%
%   Additional arguments can be supplied like this:
%   SGRAM(f,'nf','tfr',2,'log'). The arguments must be character strings
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

b=L/M;

dofs=~isempty(keyvals.fs);
Ndisp=size(coef,2);

if flags.do_nf
  % Calculate negative frequencies, use DGT
  % Display zero frequency in middle of plot.
 
  % Move zero frequency to the center.
  coef=fftshift(coef,1);

  if flags.do_fmax
    yr=-keyvals.fmax:keyvals.fmax/M:keyvals.fmax;
  else
    if dofs
      yr=-keyvals.fs/2:keyvals.fs/M:keyvals.fs/2;
    else
      yr=-L/2:b:L/2;
    end;
  end;
else
  if flags.do_fmax
    yr=0:keyvals.fmax/M:keyvals.fmax;
  else
    if dofs
      yr=0:keyvals.fs/M:keyvals.fs/2;
    else
      yr=0:b:L/2;
    end;
  end;
end;

if flags.do_tc
  xr=-floor(Ndisp/2)*a:a:floor((Ndisp-1)/2)*a;
else
  xr=0:a:Ndisp*a-1;
end;

if dofs
  % Scale x-axis by sampling rate.
  xr=xr/keyvals.fs;  
end;

% Scale x-axis by resampling rate
xr=xr/resamp;

% Apply transformation to coefficients.
if isfield(flags,'do_db') && flags.do_db
  % We have already squared the coefficients, so only multiply
  % by 10. We add realmin to avoid taking the log of zero.
  coef=10*log10(coef+realmin);
end;
  
% 'dynrange' parameter is handled by thresholding the coefficients.
if isfield(flags,'dynrange') && flags.do_dynrange
  maxclim=max(coef(:));
  coef(coef<maxclim-keyvals.dynrange)=maxclim-keyvals.dynrange;
end;

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

if flags.do_colorbar
  colorbar;
end;

axis('xy');
if ~isempty(keyvals.fs)
  xlabel('Time (s)')
  ylabel('Frequency (Hz)')
else
  xlabel('Time (samples)')
  ylabel('Frequency (samples)')
end;


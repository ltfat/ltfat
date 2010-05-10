function varargout=sgram(f,varargin)
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
  
if nargin<1
  error('Too few input arguments.');
end;

if sum(size(f)>1)>1
  error('Input must be a vector.');
end;

% Define initial value for flags and key/value pairs.
defnopos.flags.wlen={'nowlen','wlen'};
defnopos.flags.thr={'nothr','thr'};
defnopos.flags.tc={'notc','tc'};
defnopos.flags.plottype={'image','contour','mesh','pcolor'};

defnopos.flags.clim={'noclim','clim'};
defnopos.flags.fmax={'nofmax','fmax'};
defnopos.flags.log={'db','lin'};
defnopos.flags.dynrange={'nodynrange','dynrange'};
defnopos.flags.colorbar={'colorbar','nocolorbar'};

if isreal(f)
  defnopos.flags.posfreq={'posfreq','nf'};
else
  defnopos.flags.posfreq={'nf','posfreq'};
end;

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

[flags,keyvals,fs]=ltfatarghelper({'fs'},defnopos,varargin,mfilename);

% Resampling rate: Used when fmax is issued.
resamp=1;

if flags.do_tc
  f=fftshift(f);
end;

dofs=~isempty(fs);

% Downsample
if flags.do_fmax
  if dofs
    resamp=keyvals.fmax*2/fs;
  else
    resamp=keyvals.fmax*2/length(f);
  end;

  f=fftresample(f,round(length(f)*resamp));
end;

Ls=length(f);

if flags.do_posfreq
   keyvals.yres=2*keyvals.yres;
end;

[a,M,L,N,Ndisp]=gabimagepars(Ls,keyvals.xres,keyvals.yres);

% Set an explicit window length, if this was specified.
if flags.do_wlen
  keyvals.tfr=keyvals.wlen^2/L;
end;

g={'gauss',keyvals.tfr};

if flags.do_nf
  coef=abs(dgt(f,g,a,M)).^2;
else
  coef=abs(dgtreal(f,g,a,M)).^2;
end;

% Cut away zero-extension.
coef=coef(:,1:Ndisp);

if flags.do_thr
  % keep only the largest coefficients.
  coef=largestr(coef,keyvals.thr);
end

tfplot(coef,a,M,L,resamp,keyvals,flags);

if nargout>0
  varargout={coef};
end;

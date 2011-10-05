function varargout=sgram(f,varargin)
%SGRAM  Spectrogram.
%   Usage: sgram(f,op1,op2, ... );
%          sgram(f,fs,op1,op2, ... );
%          C=sgram(f, ... );
%
%   SGRAM(f) plots a spectrogram of f using a Discrete Gabor Transform (DGT).
%
%   SGRAM(f,fs) does the same for a signal with sampling rate fs (sampled
%   with fs samples per second);
%
%   SGRAM(f,fs,dynrange) additionally limits the dynamic range of the
%   plot. See the description of the 'dynrange' parameter belowe.
%
%   C=SGRAM(f, ... ) returns the image to be displayed as a matrix. Use this
%   in conjunction with IMWRITE etc. These coefficients are ONLY intended to
%   be used by post-processing image tools. Numerical Gabor signal analysis
%   and synthesis should ALWAYS be done using the DGT, IDGT, DGTREAL and
%   IDGTREAL functions.
%
%   Additional arguments can be supplied like this:
%
%C      SGRAM(f,fs,'dynrange',50).
%
%   The arguments must be character strings possibly followed by an argument:
%
%-  'dynrange',r - Limit the dynamical range to r by using a colormap in
%               the interval [chigh-r,chigh], where chigh is the highest
%               value in the plot. The default value of [] means to not
%               limit the dynamical range.
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
%-  'thr',r   - Keep only the largest fraction r of the coefficients, and
%               set the rest to zero.
%
%-  'fmax',y  - Display y as the highest frequency. Default value of []
%               means to use the Nyquest frequency.
%
%-  'xres',xres - Approximate number of pixels along x-axis /time.
%               Default value is 800
%
%-  'yres',yres - Approximate number of pixels along y-axis / frequency
%               Default value is 600
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
%   In addition to these parameteres, SGRAM accepts any of the flags
%   from NORMALIZE. The window will be normalized as specified.
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

definput.import={'ltfattranslate','normalize','tfplot'};
% Override the setting from tfplot, because SGRAM does not support the
% 'dbsq' setting (it does not make sense).
definput.flags.log={'db','lin'};

% Define initial value for flags and key/value pairs.
definput.flags.wlen={'nowlen','wlen'};
definput.flags.thr={'nothr','thr'};

if isreal(f)
  definput.flags.posfreq={'posfreq','nf'};
else
  definput.flags.posfreq={'nf','posfreq'};
end;

definput.keyvals.tfr=1;
definput.keyvals.wlen=0;
definput.keyvals.thr=0;
definput.keyvals.fmax=[];
definput.keyvals.xres=800;
definput.keyvals.yres=600;

[flags,kv,fs]=ltfatarghelper({'fs','dynrange'},definput,varargin);

% Downsample
if ~isempty(kv.fmax)
  if ~isempty(fs)
    resamp=kv.fmax*2/fs;
  else
    resamp=kv.fmax*2/length(f);
  end;

  f=fftresample(f,round(length(f)*resamp));
  kv.fs=2*kv.fmax;
end;

Ls=length(f);

if flags.do_posfreq
   kv.yres=2*kv.yres;
end;

[a,M,L,N,Ndisp]=gabimagepars(Ls,kv.xres,kv.yres);

% Set an explicit window length, if this was specified.
if flags.do_wlen
  kv.tfr=kv.wlen^2/L;
end;

g={'gauss',kv.tfr,flags.norm};

if flags.do_nf
  coef=abs(dgt(f,g,a,M));
else
  coef=abs(dgtreal(f,g,a,M));
end;

% Cut away zero-extension.
coef=coef(:,1:Ndisp);

if flags.do_thr
  % keep only the largest coefficients.
  coef=largestr(coef,kv.thr);
end

if flags.do_nf
  plotdgt(coef,a,'argimport',flags,kv);
else
  plotdgtreal(coef,a,M,'argimport',flags,kv);
end;

if nargout>0
  varargout={coef};
end;

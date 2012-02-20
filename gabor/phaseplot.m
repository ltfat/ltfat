function []=phaseplot(f,varargin)
%PHASEPLOT  Phase plot
%   Usage: phaseplot(f,op1,op2, ... );
%          phaseplot(f,fs,op1,op2, ... );
%
%   `phaseplot(f)` plots the phase of f using a |dgt|_.
%
%   `phaseplot(f,fs)` does the same for a signal with sampling rate *fs* Hz.
%
%   `phaseplot` should only be used for short signals (shorter than the
%   resolution of the screen), as there will otherwise be some visual
%   aliasing, such that very fast changing areas will look very smooth.
%   `phaseplot` always calculates the phase of the full time/frequency plane
%   (as opposed to |sgram|_), and you therefore risk running out of memory
%   for long signals.
%
%   `phaseplot` takes the following flags at the end of the line of input
%   arguments:
%
%     'tfr',v     Set the ratio of frequency resolution to time resolution.
%                 A value v=1 is the default. Setting v>1 will give better
%                 frequency resolution at the expense of a worse time
%                 resolution. A value of 0<v<1 will do the opposite.
%  
%     'wlen',s    Window length. Specifies the length of the window
%                 measured in samples. See help of PGAUSS on the exact
%                 details of the window length.
%  
%     'nf'        Display negative frequencies, with the zero-frequency
%                 centered in the middle. For real signals, this will just
%                 mirror the upper half plane. This is standard for complex
%                 signals.
%  
%     'tc'        Time centering. Move the beginning of the signal to the
%                 middle of the plot. This is usefull for visualizing the
%                 window functions of the toolbox.
%  
%     'thr',r     Keep the coefficients with a magnitude larger than r times the
%                 largest magnitude. Set the phase of the rest of the
%                 coefficients to zero. This is useful, because for small
%                 amplitude the phase values can be meaningless.
%  
%     'timeinv'   Display the phase as computed by a time-invariant
%                 |dgt|_. This is the default.
%  
%     'freqinv'   Display the phase as computed by a frequency-invariant
%                 |dgt|_.
%  
%     'fmax',y    Display y as the highest frequency.
% 
%   See also: phaselock
%
%   Demos: demo_phaseplot
%
%   References: Carmona98practical

%   AUTHOR: Peter Soendergaard
%   REFERENCE: NA
%   TESTING: NA

if nargin<1
  error('Too few input arguments.');
end;

if sum(size(f)>1)>1
  error('Input must be a vector.');
end;

% Define initial value for flags and key/value pairs.
definput.flags.wlen={'nowlen','wlen'};
definput.flags.tc={'notc','tc'};

definput.flags.clim={'noclim','clim'};
definput.flags.fmax={'nofmax','fmax'};
definput.flags.phase={'timeinv','freqinv'};

if isreal(f)
  definput.flags.posfreq={'posfreq','nf'};
else
  definput.flags.posfreq={'nf','posfreq'};
end;

definput.keyvals.fs=[];
definput.keyvals.tfr=1;
definput.keyvals.wlen=0;
definput.keyvals.thr=[];
definput.keyvals.fmax=0;

[flags,kv,fs]=ltfatarghelper({'fs'},definput,varargin);

% Resampling rate: Used when fmax is issued.
resamp=1;

if flags.do_tc
  f=fftshift(f);
end;

dofs=~isempty(fs);

% Downsample
if flags.do_fmax
  if dofs
    resamp=kv.fmax*2/fs;
  else
    resamp=kv.fmax*2/length(f);
  end;

  f=fftresample(f,round(length(f)*resamp));
end;

Ls=length(f);

% Always do the full STFT
L=length(f);
a=1;
b=1;
M=L;
N=L;

% Set an explicit window length, if this was specified.
if flags.do_wlen
  kv.tfr=kv.wlen^2/L;
end;

g={'gauss',kv.tfr};

if flags.do_nf
  coef=dgt(f,g,a,M,flags.phase);
else
  coef=dgtreal(f,g,a,M,flags.phase);
end;

if ~isempty(kv.thr)
  % keep only the largest coefficients.
  maxc=max(abs(coef(:)));
  mask=abs(coef)<maxc*kv.thr;
  coef(mask)=0;
end

coef = angle(coef);

if flags.do_nf
  % Move zero frequency to the center.
  coef=fftshift(coef,1);

  if flags.do_fmax
    yr=-fmax:fmax/M:fmax;
  else
    if dofs
      yr=-fs/2:fs/M:fs/2;
    else
      yr=-L/2:b:L/2;
    end;
  end;

else
  if flags.do_fmax
    yr=0:fmax/M:fmax;
  else
    if dofs
      yr=0:fs/M:fs/2;
    else
      yr=0:b:L/2;
    end;
  end;
end;

if flags.do_tc
  xr=-floor(N/2)*a:a:floor((N-1)/2)*a;
else
  xr=0:a:N*a-1;
end;

if ~isempty(fs)
  % Scale x-axis by sampling rate.
  xr=xr/fs;  
end;

imagesc(xr,yr,coef);
axis('xy');
if dofs
  xlabel('Time (s)')
  ylabel('Frequency (Hz)')
else
  xlabel('Time')
  ylabel('Frequency')
end;


function fftgram(f, varargin)
%FFTGRAM Plot the energy of the discrete Fourier transform
%   Usage:  fftgram(f)
%           fftgram(f, fs)
%
%   `fftgram(f)` plots the energy of the discrete Fourier transform computed 
%   from the function f. The function forms a Fourier pair with the periodic
%   autocorrelation function.
%
%   `fftgram(f,fs)` does the same for a signal sampled with a sampling
%   frequency of *fs* Hz. If *fs* is no specified, the plot will display
%   normalized frequencies.
%
%   `fftgram(f,fs,dynrange)` additionally specifies the dynamic range to
%   display on the figure.
%
%   Additional arguments for `fftgram`:
%
%      'db'      Plots the energy on a dB scale. This is the default.
%
%      'lin'     Plots the energy on a linear scale.
%
%   In addition to these parameters, `fftgram` accepts any of the flags from
%   |setnorm|. The input signal will be normalized as specified.
%
%   See also: dft, plotfft

% AUTHOR: Jordy van Velthoven

% Assert correct number of input parameters.
complainif_notenoughargs(nargin, 1, 'FFTGRAM');

definput.import={'ltfattranslate','setnorm'};
definput.keyvals.fs=[];
definput.keyvals.clim=[];
definput.keyvals.dynrange=[];  
definput.flags.powscale={'db', 'lin'};

[flags, kv] = ltfatarghelper({'fs','dynrange'},definput,varargin);

if isreal(f)
  p = (fftreal(f).*conj(fftreal(f)));
else
  p = (fft(f).*conj(fft(f)));
end;

p = setnorm(p, flags.norm);

if flags.do_db
  if isreal(f)
    plotfftreal(p,kv.fs, kv.dynrange);
  else
    plotfft(p,kv.fs, kv.dynrange);
  end;
  ylabel('Energy (dB)');
end;

if flags.do_lin
  if isreal(f)
    plotfftreal(p, kv.fs, kv.dynrange, 'lin');
  else 
    plotfft(p, kv.fs, kv.dynrange,'lin');
  end;
  ylabel('Energy');
end;

function enspecden(f, varargin)
%ENSPECDEN Plot the energy spectral density
%
%   Usage: enspecden(f);
%
%   `enspecden(f)` plots the energy spectral density of the function f.
%
%   Additional arguments for ENSPECDEN:
%
%     'db'      Plots the energy density on a dB scale.
%
%     'lin'     Plots the energy density on a linear scale.
%
%     'norm'    Normalizes the $l^2$-norm of the energy density to be 1.
%
%   See also: paucorr, dft, plotfft

% AUTHOR: Jordy van Velthoven

% Assert correct number of input parameters.
complainif_notenoughargs(nargin, 1, 'ENSPECDEN');

definput.flags.powscale={'db', 'lin'};
definput.flags.norm={'nonorm', 'norm'};

[flags, kv] = ltfatarghelper({},definput,varargin);

if isreal(f)
  p = (fftreal(f).*conj(fftreal(f)));
else
  p = (fft(f).*conj(fft(f)));
end;

if flags.do_norm
  p = normalize(p, '2');
end;


if flags.do_db
  if isreal(f)
  plotfftreal(p);
  else
  plotfft(p);
  end;
ylabel('Energy spectral density (dB)');
end;

if flags.do_lin
  if isreal(f)
  plotfftreal(p, 'lin');
  else 
  plotfft(p, 'lin');
  end;
ylabel('Energy spectral density');
end;
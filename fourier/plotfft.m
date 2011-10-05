function plotfft(coef,varargin)
%PLOTFFT  Plot the output from FFT
%   Usage: plotfft(coef);
%          plotfft(coef,fs);
%
%   PLOTFFT(coef) will plot the output from the FFT function. The
%   frequency axis will use normalized frequencies between 0 and 1 (Nyquest).
%
%   PLOTFFT(coef,fs) will do the same for the FFT of a signal sampled at
%   a sampling rate of fs Hz.
%
%   See also: plotfftreal
  
if nargin<1
  error('%s: Too few input parameters.',upper(mfilename));
end;

if ~isvector(coef)>1
  error('Input is multidimensional.');
end;

definput.import={'ltfattranslate','normalize'};
definput.importdefaults={'null'};

definput.flags.log={'db','dbsq','lin','linsq','linabs'};
definput.flags.posfreq={'nf','posfreq'};

definput.keyvals.fs=[];
definput.keyvals.dynrange=[];

definput.keyvals.opts={};

[flags,kv,fs]=ltfatarghelper({'fs','dynrange'},definput,varargin);

coef=normalize(coef,flags.norm);

if ~isvector(coef)
  error('%s: Can only plot vectors.',upper(mfilename));
end;
N=length(coef);

% Apply transformation to coefficients.
if flags.do_db
  coef=20*log10(abs(coef)+realmin);
end;

if flags.do_dbsq
  coef=10*log10(abs(coef)+realmin);
end;

if flags.do_linsq
  coef=abs(coef).^2;
end;

if flags.do_linabs
  coef=abs(coef);
end;

if flags.do_lin
  if ~isreal(coef)
    error(['Complex valued input cannot be plotted using the "lin" flag.',...
           'Please use the "linsq" or "linabs" flag.']);
  end;
end;
  
% 'dynrange' parameter is handled by thresholding the coefficients.
if ~isempty(kv.dynrange)
  maxclim=max(coef(:));
  coef(coef<maxclim-kv.dynrange)=maxclim-kv.dynrange;
end;

if flags.do_nf
  if rem(N,2)==0
    xr=(-N/2+1:N/2)*2/N;  
    % Subtract 1 in order to place the Nyquest frequency following the
    % positive frequencies. That is why we do not use fftshift.
    coef=circshift(coef,N/2-1);
  else
    xr=(-floor(N/2):floor(N/2))*2/N;  
    coef=fftshift(coef);
  end;
else

  N2=floor(N/2)+1;
  coef=coef(1:N2);
  xr=(0:N2-1)*2/N;  
end;

if ~isempty(kv.fs)
  xr=xr*kv.fs/2;
end;

plot(xr,coef,kv.opts{:});
xlim([xr(1) xr(end)]);

ylabel(sprintf('%s (dB)',kv.magnitude));
if ~isempty(kv.fs)
  xlabel(sprintf('%s (Hz)',kv.frequency));
else
  xlabel(sprintf('%s (%s)',kv.frequency,kv.normalized));
end;

function plotfftreal(coef,varargin)  
%PLOTFFTREAL  Plot the output from FFTREAL  
%   Usage: plotfftreal(coef);
%          plotfftreal(coef,N);
%          plotfftreal(coef,N,fs);
%
%   `plotfftreal(coef)` plots the output from the |fftreal|_ function. The
%   frequency axis will use normalized frequencies between 0 and 1 (the
%   Nyquest frequency). It is assumed that the length of the original
%   transform was even.
%
%   `plotfftreal(coef,N)` additionally specifies then length of the
%   original transform. Use this if you are unsure if the transform
%   length was even.
%
%   `plotfft(coef,N,fs)` does the same for the |fftreal|_ of a signal
%   sampled at a sampling rate of *fs* Hz.
%
%   `plotfftreal` accepts the following optional arguments:
%
%     'db'     Apply $20\cdot \log_{10}$ to the coefficients. This makes 
%              it possible to see very weak phenomena, but it might show 
%              too much noise. This is the default.
%
%     'dbsq'   Apply $10\cdot \log_{10}$ to the coefficients. Same as the
%              `'db'` option, but assumes that the input is already squared.  
%
%     'lin'    Show the coefficients on a linear scale. This will
%              display the raw input without any modifications. Only works for
%              real-valued input.
%
%     'linsq'  Show the square of the coefficients on a linear scale.
%
%     'linabs'  Show the absolute value of the coefficients on a linear scale.
%
%   In addition to these parameteres, `plotfftreal` accepts any of the flags
%   from |normalize|_. The coefficients will be normalized as specified
%   before plotting.
%
%   See also: plotfft, fftreal

  
if nargin<1
  error('%s: Too few input parameters.',upper(mfilename));
end;

if ~isvector(coef)>1
  error('Input is multidimensional.');
end;

definput.import={'ltfattranslate','normalize'};
definput.importdefaults={'null'};

definput.flags.log={'db','dbsq','lin','linsq','linabs'};

definput.keyvals.fs=[];
definput.keyvals.dynrange=[];
definput.keyvals.opts={};

definput.keyvals.N=2*(length(coef)-1);
[flags,kv,N,fs]=ltfatarghelper({'N','fs','dynrange'},definput,varargin);

N2=floor(N/2)+1;
if N2~=length(coef)
  error('%s: Size mismatch.',upper(mfilename));
end;

coef=normalize(coef,flags.norm);

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

xr=(0:N2-1)*2/N;
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


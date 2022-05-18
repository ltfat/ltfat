function plotfftreal(coef,varargin)  
%PLOTFFTREAL  Plot the output from FFTREAL  
%   Usage: plotfftreal(coef);
%          plotfftreal(coef,fs);
%
%   `plotfftreal(coef)` plots the output from the |fftreal| function. The
%   frequency axis will use normalized frequencies between 0 and 1 (the
%   Nyquist frequency). It is assumed that the length of the original
%   transform was even.
%
%   `plotfftreal(coef,fs)` does the same for the |fftreal| of a signal
%   sampled at a sampling rate of *fs* Hz.
%
%   `plotfftreal(coef,fs,dynrange)` additionally limits the dynamic range of the
%   plot. See the description of the `'dynrange'` parameter below.
%
%   `plotfftreal` accepts the following optional arguments:
%
%     'dynrange',r  Limit the dynamical range to `r` by using a colormap in
%                   the interval `[chigh-r,chigh]`, where `chigh` is the highest
%                   value in the plot. The default value of `[]` means to not
%                   limit the dynamical range. 
%
%     'db'      Apply $20\cdot \log_{10}$ to the coefficients. This makes 
%               it possible to see very weak phenomena, but it might show 
%               too much noise. This is the default.
%
%     'dbsq'    Apply $10\cdot \log_{10}$ to the coefficients. Same as the
%               `'db'` option, but assumes that the input is already squared.  
%
%     'lin'     Show the coefficients on a linear scale. This will
%               display the raw input without any modifications. Only works for
%               real-valued input.
%
%     'linsq'   Show the square of the coefficients on a linear scale.
%
%     'linabs'  Show the absolute value of the coefficients on a linear
%               scale.
%     
%     'N',N     Specify the transform length *N*. Use this if you are
%               unsure if the original input signal was of even length.
%
%     'dim',dim  If `coef` is multidimensional, `dim` indicates the 
%                dimension along which are the individual channels oriented.
%                Value 1 indicates columns, value 2 rows.
%
%     'flog'  Use logarithmic scale for the frequency axis.
%
%
%   In addition to these parameters, `plotfftreal` accepts any of the flags
%   from |setnorm|. The coefficients will be normalized as specified
%   before plotting.
%
%   See also: plotfft, fftreal

  
if nargin<1
  error('%s: Too few input parameters.',upper(mfilename));
end;

definput.import={'ltfattranslate','setnorm'};
definput.importdefaults={'null'};

definput.flags.log={'db','dbsq','lin','linsq','linabs'};
definput.flags.freqscale={'flin','flog'};


definput.keyvals.fs=[];
definput.keyvals.dynrange=[];
definput.keyvals.opts={};

definput.keyvals.N=[];
definput.keyvals.dim=[];
[flags,kv,fs]=ltfatarghelper({'fs','dynrange'},definput,varargin);

% if ~isvector(coef)
%   error('%s: Input is multidimensional.',upper(mfilename));
% end;

[coef,~,Lc]=assert_sigreshape_pre(coef,[],kv.dim,upper(mfilename));

N=kv.N;
if isempty(N)
   N=2*(Lc-1);
end

N2=floor(N/2)+1;
if N2~=Lc
  error('%s: Size mismatch.',upper(mfilename));
end;

coef=setnorm(coef,flags.norm);

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

if flags.do_flin
   plot(xr,coef,kv.opts{:});
elseif flags.do_flog
   semilogx(xr,coef,kv.opts{:}); 
end
xlim([xr(1) xr(end)]);


if flags.do_db || flags.do_dbsq
  ylabel(sprintf('%s (dB)',kv.magnitude));
else
  ylabel(sprintf('%s',kv.magnitude));
end;

if ~isempty(kv.fs)
  xlabel(sprintf('%s (Hz)',kv.frequency));
else
  xlabel(sprintf('%s (%s)',kv.frequency,kv.normalized));
end;


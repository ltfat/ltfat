function plotfftreal(coef,varargin)  
%PLOTFFTREAL  Plot the output from FFTREAL  
%   Usage: plotfftreal(coef);
%          plotfftreal(coef,N);
%          plotfftreal(coef,
if nargin<1
  error('%s: Too few input parameters.',upper(mfilename));
end;

if ~isvector(coef)>1
  error('Input is multidimensional.');
end;

definput.import={'normalize'};
definput.importdefaults={'null'};

definput.flags.log={'db','dbsq','lin','linsq','linabs'};

definput.keyvals.fs=[];
definput.keyvals.dynrange=[];

definput.keyvals.frequency='Frequency';
definput.keyvals.time='Time';
definput.keyvals.samples='samples';
definput.keyvals.normalized='normalized';

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
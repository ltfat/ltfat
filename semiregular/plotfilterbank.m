function [] = plotfilterbank(coef,a,varargin)
%PLOTNSDGT Plot spectrogram from filterbank coefficients
%   Usage:  plotfilterbank(coef,a);
%           plotfilterbank(coef,a,fc);
%           plotfilterbank(coef,a,fc,fs);
%           plotfilterbank(coef,a,fc,fs,dynrange);
%
%   PLOTFILTERBANK(coef,a) plots filterbank coefficients coef obtained from
%   either the FILTERBANK or UFILTERBANK functions. The coefficients must
%   have been produced with a time-shift of _a. For more details on the
%   format of the variables coef and _a, see the help of the FILTERBANK
%   or UFILTERBANK functions.
%
%   PLOTFILTERBANK(coef,a,fc) makes it possible to specify the center
%   frequency for each channel in the vector fc.
%
%   PLOTFILTERBANK(coef,a,fc,fs) does the same assuming a sampling rate of
%   fs Hz of the original signal.
%
%   PLOTFILTERBANK(coef,a,fc,fs,dynrange) makes it possible to specify
%   the dynamic range of the coefficients.
%
%   PLOTFILTERBANK supports all the optional parameters of TFPLOT. Please
%   see the help of TFPLOT for an exhaustive list.
%
%   In addition to the flags and key/values in TFPLoT, PLOTFILTERBANK
%   supports the following optional arguments:
%
%-     'fc',fc    - Center frequencies of the channels. fc must be a vector with
%                   the length equal to the number of channels. The
%                   default value of [] means to plot the channel
%                   no. instead of its frequency.
%
%-     'ntickpos',n - Number of tick positions along the y-axis. The
%                   position of the ticks are determined automatically.
%                   Default value is 10.
%
%-     'tick',t     - Array of tick positions on the y-axis. Use this
%                   option to specify the tick position manually.
%
%   See also:  filterbank, ufilterbank, tfplot, sgram

if nargin<2
  error('%s: Too few input parameters.',upper(mfilename));
end;

definput.import={'tfplot'};
definput.keyvals.fc=[];
definput.keyvals.ntickpos=10;
definput.keyvals.tick=[];

[flags,kv,fs]=ltfatarghelper({'fc','fs','dynrange'},definput,varargin);
  
N=size(coef,1);
M=size(coef,2);

% Turn the coefficients as in DGT.
coef=coef.';

% Freq. pos is just number of the channel.
yr=1:M;

if size(coef,3)>1
  error('Input is multidimensional.');
end;

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
  
% 'dynrange' parameter is handled by thresholding the coefficients.
if ~isempty(kv.dynrange)
  maxclim=max(coef(:));
  coef(coef<maxclim-kv.dynrange)=maxclim-kv.dynrange;
end;

if flags.do_tc
  xr=(-floor(N/2):floor((N-1)/2))*a;
else
  xr=(0:N-1)*a;
end;

if ~isempty(kv.fs)
  xr=xr/kv.fs;
end;

switch(flags.plottype)
  case 'image'
    if flags.do_clim
      imagesc(xr,yr,coef,kv.clim);
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
if ~isempty(kv.fs)
  xlabel(sprintf('%s (s)',kv.time));
else
  xlabel(sprintf('%s (%s)',kv.time,kv.samples));
end;

if isempty(kv.fc)
  ylabel('Channel No.');
else
  
  if isempty(kv.tick)
    tickpos=linspace(1,M,kv.ntickpos);
    tick=spline(1:M,kv.fc,tickpos);
    
    set(gca,'YTick',tickpos);
    set(gca,'YTickLabel',num2str(tick(:),3));

  else
    % Create a crude inverse mapping to determine the positions of the
    % ticks
    tick=kv.tick;
    manyticks=spline(1:M,kv.fc,linspace(1,M,1000));
    nticks=length(tick);
    tickpos=zeros(nticks,1);
    for ii=1:nticks
      jj=find(manyticks>=tick(ii));
      tickpos(ii)=jj(1)/1000*M;
    end;
    
    set(gca,'YTick',tickpos);
    set(gca,'YTickLabel',num2str(tick(:)));

  end;
  
  ylabel(sprintf('%s (s)',kv.frequency));


  
end;


function [] = plotfilterbank(coef,a,varargin)
%PLOTFILTERBANK Plot filterbank and ufilterbank coefficients
%   Usage:  plotfilterbank(coef,a);
%           plotfilterbank(coef,a,fc);
%           plotfilterbank(coef,a,fc,fs);
%           plotfilterbank(coef,a,fc,fs,dynrange);
%
%   `plotfilterbank(coef,a)` plots filterbank coefficients *coef* obtained from
%   either the |filterbank|_ or |ufilterbank|_ functions. The coefficients must
%   have been produced with a time-shift of *a*. For more details on the
%   format of the variables *coef* and *a*, see the help of the |filterbank|_
%   or |ufilterbank|_ functions.
%
%   `plotfilterbank(coef,a,fc)` makes it possible to specify the center
%   frequency for each channel in the vector *fc*.
%
%   `plotfilterbank(coef,a,fc,fs)` does the same assuming a sampling rate of
%   *fs* Hz of the original signal.
%
%   `plotfilterbank(coef,a,fc,fs,dynrange)` makes it possible to specify
%   the dynamic range of the coefficients.
%
%   `plotfilterbank` supports all the optional parameters of |tfplot|_. Please
%   see the help of |tfplot|_ for an exhaustive list.
%
%   In addition to the flags and key/values in |tfplot|_, `plotfilterbank`
%   supports the following optional arguments:
%
%     'fc',fc       Centre frequencies of the channels. *fc* must be a vector with
%                   the length equal to the number of channels. The
%                   default value of [] means to plot the channel
%                   no. instead of its frequency.
%
%     'ntickpos',n  Number of tick positions along the y-axis. The
%                   position of the ticks are determined automatically.
%                   Default value is 10.
%
%     'tick',t      Array of tick positions on the y-axis. Use this
%                   option to specify the tick position manually.
%
%     'audtick'     Use ticks suitable for visualizing an auditory
%                   filterbank. Same as `'tick',[0,100,250,500,1000,...]`.
%
%   See also:  filterbank, ufilterbank, tfplot, sgram

if nargin<2
  error('%s: Too few input parameters.',upper(mfilename));
end;

definput.import={'plotfilterbank','tfplot','ltfattranslate'};

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
    nlarge=1000;
    tick=kv.tick;
        
    % Create a crude inverse mapping to determine the positions of the
    % ticks. Include half a channel in each direction, because it is
    % possible to display a tick mark all the way to the edge of the
    % plot.
    manyticks=spline(1:M,kv.fc,linspace(0.5,M+0.5,nlarge));
    
    % Keep only ticks <= than highest frequency+.5*bandwidth
    tick=tick(tick<=manyticks(end));
    
    % Keep only ticks >= lowest frequency-.5*bandwidth
    tick=tick(tick>=manyticks(1));    
    
    nticks=length(tick);
    tickpos=zeros(nticks,1);
    for ii=1:nticks
      jj=find(manyticks>=tick(ii));
      tickpos(ii)=jj(1)/nlarge*M;
    end;

    set(gca,'YTick',tickpos);
    set(gca,'YTickLabel',num2str(tick(:)));

  end;
  
  ylabel(sprintf('%s (Hz)',kv.frequency));


  
end;


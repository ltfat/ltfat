function [ ] = plotfwt(c,varargin )
%PLOTFWT  Plot wavelet coefficients
%   Usage:  plotfwt(c) 
%           plotfwt(c,Lc)
%           plotfwt(...,type,fs,dynrange)
%
%   `plotfwt(c)` plots the wavelet coefficients *c*. The input cell-array
%   *c* must have the same structure as coefficients returned from
%   |fwt|_. See the help on |fwt|_ for more information.
%
%   `plotfwt(c,fs)` does the same assuming a sampling rate of *fs* Hz of the
%   original signal.
%
%   `plotfwt(c,fs,dynrange)` additionally limits the dynamic range.
%
%   `C=plotfwt(...)` returns the processed image data used in the
%   plotting. Inputting this data directly to `imagesc` or similar functions
%   will create the plot. This is usefull for custom post-processing of the
%   image data.
%
%   `plotfwt` supports all the optional parameters of |tfplot|_. Please see
%   the help of |tfplot|_ for an exhaustive list.
%
%   See also: fwt, tfplot

if nargin<1
  error('%s: Too few input parameters.',upper(mfilename));
end;

definput.import={'tfplot','ltfattranslate','fwt','plotfwt'};

definput.keyvals.xres=800;

[flags,kv]=ltfatarghelper({'fs','dynrange'},definput,varargin);


if iscell(c)
  cM = length(c); 
  L=size(c{end}(:,1),1); 


  N=kv.xres;
  if(flags.do_uni)
     coef2=zeros(cM,N);
     yr=1:cM;
  elseif(flags.do_dyad)
     coef2=zeros(2^(cM-1),N); 
     yr=1:2^(cM-1);
  end
  
  row=c{1}(:,1);
  coef2(end,:)=interp1(linspace(0,1,numel(row)),row,linspace(0,1,N),'nearest');
  runii = 1;
  for ii=1:cM-1
       row=c{end+1-ii}(:,1);
       if(flags.do_uni)
            rows = 1;
       elseif(flags.do_dyad)
            rows = 2^(cM-ii-1);
       end

    rowint = interp1(linspace(0,1,numel(row)),row,linspace(0,1,N),'nearest');

    coef2(runii:runii+rows-1,:)= repmat(rowint,[rows,1]);
    runii = runii + rows;
  end;
  coef=coef2;
  delta_t=L/N;
else
 
end;

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
 
% Handle clim by thresholding and cutting
if ~isempty(kv.clim)
  coef(coef<kv.clim(1))=kv.clim(1);
  coef(coef>kv.clim(2))=kv.clim(2);
end;
 


if flags.do_tc
  xr=(-floor(N/2):floor((N-1)/2))*a;
else
  xr=(0:N-1)*delta_t;
end;

if ~isempty(kv.fs)
  xr=xr/kv.fs;
end;

if flags.do_display

    switch(flags.wavplottype)
      case 'image'
        imagesc(xr,yr,coef);
        hold on;
        
        hold off;
      case 'waterfall'
        waterfall(xr,yr,coef);
        hold on;
        
        hold off; 
      case 'stem'
        subplot(M,1,1);
        stem(c{1});
        axis tight;
        for jj=2:cM
            for ff=1:cN
                subplot(M,1,1+(jj-2)*cN+ff);
                stem(c{jj,ff});
                axis tight;
            end;
        end;
      case 'surf'
        surf(xr,yr,coef,'EdgeColor','none');
    end;
    
    if flags.do_colorbar
        colorbar;
    end;
    
    if ~isempty(kv.fs)
        xlabel(sprintf('%s (s)',kv.time));
        ylabel(sprintf('%s (Hz)',kv.frequency));
    else
        xlabel(sprintf('%s (%s)',kv.time,kv.samples));
        ylabel('Level XXX');
    end;

end;
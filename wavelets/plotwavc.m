function [ ] = plotwavc(c,varargin )
%PLOTWAVC  Plot wavelet coefficients
%   Usage:  plotwavc(c) 
%           plotwavc(c,Lc)
%           plotwavc(...,type,fs,dynrange)
%
%   `plotwavc(c)` plot the wavelet coefficients *c*. The input cell-array
%   *c* must have the same structure as coefficients returned from
%   |fwt|_. See the help on |fwt|_ for more information.
%
%   See also: fwt

if nargin<1
  error('%s: Too few input parameters.',upper(mfilename));
end;

definput.import={'tfplot','ltfattranslate','fwt','plotwavc'};

definput.keyvals.xres=800;

[flags,kv]=ltfatarghelper({'fs','dynrange'},definput,varargin);





if iscell(c)
  [cM, cN] = size(c); 

  M=numel(c);
  if cN>1
      M=M-(cN-1);
  end

     if(flags.do_undec)
       L=size(c{1},1); 
     else
       L=size(c{1},1)*2^(cM-1)+2^(cM-1);  
     end

  N=kv.xres;
  if(flags.do_uni)
     coef2=zeros(M,N);
     yr=1:M;
  elseif(flags.do_dyad)
     coef2=zeros(2^(M-1),N); 
     yr=1:2^(M-1);
  end
  
  row=c{1}(:,1);
  coef2(end,:)=interp1(linspace(0,1,numel(row)),row,linspace(0,1,N),'nearest');
  runii = 1;
  for ii=1:cM-1
    for ff=1:cN
       row=c{end+1-ii,end+1-ff}(:,1);
       if(flags.do_uni)
            rows = 1;
       elseif(flags.do_dyad)
            rows = 2^(M-ii-1);
       end

    rowint = interp1(linspace(0,1,numel(row)),row,linspace(0,1,N),'nearest');

    coef2(runii:runii+rows-1,:)= repmat(rowint,[rows,1]);
    runii = runii + rows;
    end;
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
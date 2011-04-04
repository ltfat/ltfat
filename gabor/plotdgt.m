function plotdgt(coef,a,varargin)
%DGTPLOT  Plot DGT coefficients.
%   Usage: plotdgt(coef,a);
%          plotdgt(coef,a,fs);
%          plotdgt(coef,a,fs,dynrange);
%
%   See also:  sgram

%   AUTHOR : Peter Soendergaard.
%   TESTING: NA
%   REFERENCE: NA

if nargin<2
  error('%s: Too few input parameters.',upper(mfilename));
end;

definput.flags.wlen={'nowlen','wlen'};
definput.flags.thr={'nothr','thr'};
definput.flags.tc={'notc','tc'};
definput.flags.plottype={'image','contour','mesh','pcolor'};
definput.flags.clim={'noclim','clim'};
definput.flags.log={'db','lin'};
definput.flags.envtype={'abs','sq'};
definput.flags.colorbar={'colorbar','nocolorbar'};

definput.keyvals.fs=[];
definput.keyvals.thr=0;
definput.keyvals.clim=[0,1];
definput.keyvals.dynrange=[];
definput.keyvals.xres=800;

[flags,kv,fs]=ltfatarghelper({'fs','dynrange'},definput,varargin);


M=size(coef,1);
N=size(coef,2);
L=N*a;
b=L/M;

if ~isreal(coef)
  coef=abs(coef);
end;

dofs=~isempty(kv.fs);
Ndisp=size(coef,2);

% Move zero frequency to the center.
coef=fftshift(coef,1);

if dofs
  yr=-kv.fs/2:kv.fs/M:kv.fs/2;
else
  %yr=-L/2:b:L/2;      
  yr=linspace(-1,1,M);
end;

if flags.do_tc
  xr=-floor(Ndisp/2)*a:a:floor((Ndisp-1)/2)*a;
else
  xr=0:a:Ndisp*a-1;
end;

if dofs
  % Scale x-axis by sampling rate.
  xr=xr/kv.fs;  
end;

% Apply transformation to coefficients.
if flags.do_db
  % We have already squared the coefficients, so only multiply
  % by 10. We add realmin to avoid taking the log of zero.
  if flags.do_sq
    coef=10*log10(coef+realmin);
  else
    coef=20*log10(coef+realmin);
  end;
end;
  
% 'dynrange' parameter is handled by thresholding the coefficients.
if ~isempty(kv.dynrange)
  maxclim=max(coef(:));
  coef(coef<maxclim-kv.dynrange)=maxclim-kv.dynrange;
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
  xlabel('Time (s)')
  ylabel('Frequency (Hz)')
else
  xlabel('Time (samples)')
  ylabel('Frequency (normalized)')
end;


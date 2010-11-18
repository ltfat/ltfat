function tfplot(coef,a,M,L,resamp,keyvals,flags)
%TFPLOT  Healper routine for SGRAM and other TF plotting routines.
%   Usage: tfplot(coef,a,M,L,resamp,keyvals,flags)
%
%   Don't call this function directly, see the code for SGRAM on how to
%   use it.
%
%   See also:  sgram

%   AUTHOR : Peter Soendergaard.
%   TESTING: NA
%   REFERENCE: NA

b=L/M;

dofs=~isempty(keyvals.fs);
Ndisp=size(coef,2);

if flags.do_nf
  % Calculate negative frequencies, use DGT
  % Display zero frequency in middle of plot.
 
  % Move zero frequency to the center.
  coef=fftshift(coef,1);

  if flags.do_fmax
    yr=-keyvals.fmax:keyvals.fmax/M:keyvals.fmax;
  else
    if dofs
      yr=-keyvals.fs/2:keyvals.fs/M:keyvals.fs/2;
    else
%yr=-L/2:b:L/2;      
      yr=linspce(-1,1,M);
    end;
  end;
else
  if flags.do_fmax
    yr=0:keyvals.fmax/M:keyvals.fmax;
  else
    if dofs
      yr=0:keyvals.fs/M:keyvals.fs/2;
    else
      %yr=0:b:L/2;
      yr=linspace(0,1,M/2);
    end;
  end;
end;

if flags.do_tc
  xr=-floor(Ndisp/2)*a:a:floor((Ndisp-1)/2)*a;
else
  xr=0:a:Ndisp*a-1;
end;

if dofs
  % Scale x-axis by sampling rate.
  xr=xr/keyvals.fs;  
end;

% Scale x-axis by resampling rate
xr=xr/resamp;

% Apply transformation to coefficients.
if isfield(flags,'do_db') && flags.do_db
  % We have already squared the coefficients, so only multiply
  % by 10. We add realmin to avoid taking the log of zero.
  coef=10*log10(coef+realmin);
end;
  
% 'dynrange' parameter is handled by thresholding the coefficients.
if isfield(flags,'dynrange') && flags.do_dynrange
  maxclim=max(coef(:));
  coef(coef<maxclim-keyvals.dynrange)=maxclim-keyvals.dynrange;
end;

switch(flags.plottype)
  case 'image'
    if flags.do_clim
      imagesc(xr,yr,coef,keyvals.clim);
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
if ~isempty(keyvals.fs)
  xlabel('Time (s)')
  ylabel('Frequency (Hz)')
else
  xlabel('Time (samples)')
  ylabel('Frequency (normalized)')
end;


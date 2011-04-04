function tfplot(coef,step,yr,varargin)
%DGTPLOT  Plot DGT coefficients.
%   Usage: plotdgt(coef,a);
%          plotdgt(coef,a,fs);
%          plotdgt(coef,a,fs,dynrange);
%
%   PLOTDGT(coef,a) will plot the Gabor coefficient coefficients
%   coef. The coefficient must have been produce with a timeshift of _a.
%
%   The arguments must be character strings possibly followed by an argument:
%
%-  'dynrange',r - Limit the dynamical range to r by using a colormap in
%               the interval [chigh-r,chigh], where chigh is the highest
%               value in the plot. The default value of [] means to not
%               limit the dynamical range.
%
%-  'db'      - Apply 20*log10 to the coefficients. This makes it possible to
%               see very weak phenomena, but it might show too much noise. A
%               logarithmic scale is more adapted to perception of sound.
%               This is the default.
%
%-  'dbsq'    - Apply 10*log10 to the coefficients. Same as the 'db'
%               option, but assume that the input is already squared.  
%
%-  'lin'     - Show the coefficients on a linear scale.
%
%-  'tc'      - Time centering. Move the beginning of the signal to the
%               middle of the plot. 
%
%-  'clim',[clow,chigh] - Use a colormap ranging from clow to chigh. These
%               values are passed to IMAGESC. See the help on IMAGESC.
%
%-  'image'   - Use 'imagesc' to display the plot. This is the default.
%
%-  'contour' - Do a contour plot.
%          
%-  'surf'    - Do a surf plot.
%
%-  'colorbar' - Display the colorbar. This is the default.
%
%-  'nocolorbar' - Do not display the colorbar.
%
%   See also:  dgtrealplot, sgram

%   AUTHOR : Peter Soendergaard.
%   TESTING: NA
%   REFERENCE: NA

if nargin<2
  error('%s: Too few input parameters.',upper(mfilename));
end;

definput.import={'tfplot'};
[flags,kv,fs]=ltfatarghelper({'fs','dynrange'},definput,varargin);

M=size(coef,1);
N=size(coef,2);

if ~isreal(coef)
  coef=abs(coef);
end;

if size(coef,3)>1
  error('Input is multidimensional.');
end;

% Apply transformation to coefficients.
if flags.do_db
  coef=20*log10(coef+realmin);
end;

if flags.do_dbsq
  coef=10*log10(coef+realmin);
end;
  
% 'dynrange' parameter is handled by thresholding the coefficients.
if ~isempty(kv.dynrange)
  maxclim=max(coef(:));
  coef(coef<maxclim-kv.dynrange)=maxclim-kv.dynrange;
end;

if flags.do_tc
  xr=-floor(N/2)*step:step:floor((N-1)/2)*step;
else
  xr=0:step:N*step-1;
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


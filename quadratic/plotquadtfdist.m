function plotquadtfdist(p, varargin);
%PLOTQUADTFDIST Plot quadratic time-frequency distribution
%   Usage: plotquadtfdist(p);
% 
%   'plotquadtfdist(p)' plots the quadratic time-frequency distribution 
%   on the time-frequency plane. The quadratic time-frequency distribution
%   should be a square matrix.
%
%   PLOTQUADTFDIST takes the following additional arguments:
%
%     'dynrange',r
%              Limit the dynamical range to r by using a colormap in
%              the interval [chigh-r,chigh], where chigh is the highest
%              value in the plot. The default value of [] means to not
%              limit the dynamical range.
%
%     'db'     Apply 20*log_{10} to the coefficients. This makes 
%              it possible to see very weak phenomena, but it might show 
%              too much noise. A logarithmic scale is more adapted to 
%              perception of sound. This is the default.
%
%     'dbsq'   Apply 10*log_{10} to the coefficients. Same as the
%              'db' option, but assume that the input is already squared.  
%
%     'lin'    Show the coefficients on a linear scale. This will
%              display the raw input without any modifications. Only works for
%              real-valued input.
%
%     'linsq'  Show the square of the coefficients on a linear scale.
%
%     'linabs'  Show the absolute value of the coefficients on a linear scale.
%
%     'tc'     Time centring. Move the beginning of the signal to the
%              middle of the plot. 
%
%     'clim',clim   Use a colormap ranging from clim(1) to clim(2). These
%                   values are passed to imagesc. See the help on imagesc.
%
%     'image'       Use imagesc to display the plot. This is the default.
%
%     'contour'     Do a contour plot.
%          
%     'surf'        Do a surf plot.
%
%     'colorbar'    Display the colorbar. This is the default.
%
%     'nocolorbar'  Do not display the colorbar.
%
%     'display'     Display the figure. This is the default.
%

% AUTHOR: Jordy van Velthoven

complainif_notenoughargs(nargin, 1, 'PLOTQUADTFDIST');

[N1, N2] = size(p);

if N1==N2
 yr = [0, 1-2/N1];
else
 error('%s: The input should be a square matrix.', upper(mfilename));
end;

tfplot(p, 1, yr, varargin{:});


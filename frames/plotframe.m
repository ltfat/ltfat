function outsig=plotframe(insig,F,varargin);
%PLOTFRAME  Plot frame coefficients
%   Usage: plotframe(c,f);
%
%   `plotframe(c,F)` plots the frame coefficients *c* using the plot
%   command associanted to the frame *F*.
%
%   `plotframe(c,F,...)` passes any additional parameters to the native
%   plot routine. Please see the help on the specific plot routine for a
%   complete description. 
%
%   The following common set of parameters are supported by all plotting
%   routines:
%
%     'dynrange',r
%              Limit the dynamical range to *r*. The default value of []
%              means to not limit the dynamical range.
%
%     'db'     Apply $20\cdot \log_{10}$ to the coefficients. This makes 
%              it possible to see very weak phenomena, but it might show 
%              too much noise. A logarithmic scale is more adapted to 
%              perception of sound. This is the default.
%
%     'dbsq'   Apply $10\cdot \log_{10}$ to the coefficients. Same as the
%              `'db'` option, but assume that the input is already squared.  
%
%     'lin'    Show the coefficients on a linear scale. This will
%              display the raw input without any modifications. Only works for
%              real-valued input.
%
%     'linsq'  Show the square of the coefficients on a linear scale.
%
%     'linabs'  Show the absolute value of the coefficients on a linear scale.
%
%     'clim',clim
%              Only show values in between $clim(1)$ and $clim(2)$. This
%              is usually done by adjusting the colormap. See the help on `imagesc`.
%
%   See also: newframe, framet
  
switch(F.type)
 case 'dgt'
  [MN,W]=size(insig);
  N=MN/F.M;
  insig=reshape(insig,[F.M,N,W]);

  plotdgt(insig,F.a,varargin);
 case 'dgtreal'
  [MN,W]=size(insig);
  M2=floor(F.M/2)+1;
  N=MN/M2;
  insig=reshape(insig,[M2,N,W]);
  plotdgtreal(insig,F.a,F.M,varargin);
 case 'dwilt'
  [MN,W]=size(insig);
  N=MN/F.M;
  insig=reshape(insig,[2*F.M,N/2,W]);
  plotdwilt(insig,varargin);
 case 'wmdct'
  [MN,W]=size(insig);
  N=MN/F.M;
  insig=reshape(insig,[2*F.M,N/2,W]);
  plotwmdct(insig,varargin);
end;

  
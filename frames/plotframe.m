function outsig=plotframe(F,insig,varargin);
%PLOTFRAME  Plot frame coefficients
%   Usage: plotframe(F,insig,...);
%
%   `plotframe(F,c)` plots the frame coefficients *c* using the plot
%   command associated to the frame *F*.
%
%   `plotframe(F,c,...)` passes any additional parameters to the native
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
%   See also: frame, frana

switch(F.type)
   case {'fwt','ufwt','wfbt','wpfbt','uwfbt','uwpfbt'}
      info.fname = F.type;
      info.wt = F.g;
end;

switch(F.type)
 case 'dgt'
  plotdgt(framecoef2native(F,insig),F.a,varargin{:}); 
 case 'dgtreal'
  outsig = plotdgtreal(framecoef2native(F,insig),F.a,F.M,varargin{:}); 
 case 'dwilt'
  plotdwilt(framecoef2native(F,insig),varargin{:}); 
 case 'wmdct'
  plotwmdct(framecoef2native(F,insig),varargin{:});
 case 'gen'
  error(['%s: There is no default way of visualizing general frame ' ...
         'coefficients.'],upper(mfilename));
 case 'dft'
  plotfft(insig,varargin{:});
 case {'dcti','dctii','dctiii','dctiv',...
       'dsti','dstii','dstiii','dstiv'}
  % FIXME : This is not strictly correct, as half the transforms use an
  % odd frequency centering.
  plotfftreal(insig,varargin{:});
 case 'fwt'
    info.Lc = fwtclength(size(insig,1)/F.red,F.g,F.J);
    info.J = F.J;
    info.dim = 1;
    plotwavelets(insig,info,varargin{:});  
 case 'ufwt'
    info.J = F.J;
    plotwavelets(framecoef2native(F,insig),info,varargin{:}); 
 case {'wfbt','wpfbt'}
    plotwavelets(framecoef2native(F,insig),info,varargin{:}); 
 case {'uwfbt','uwpfbt'}
    plotwavelets(framecoef2native(F,insig),info,varargin{:}); 
    
case {'filterbank','filterbankreal'}
    plotfilterbank(framecoef2native(F,insig),F.a,[],varargin{:});
end;

  
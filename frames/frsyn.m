function outsig=frsyn(F,insig);
%FRSYN  Frame synthesis operator
%   Usage: f=frsyn(F,c);
%
%   `f=frsyn(F,c)` constructs a signal *f* from the frame coefficients *c*
%   using the frame *F*. The frame object *F* must have been created using
%   |newframe|_.
%
%   See also: newframe, frana, plotframe
  
if nargin<2
  error('%s: Too few input parameters.',upper(mfilename));
end;

if ~isstruct(F)
  error('%s: First agument must be a frame definition structure.',upper(mfilename));
end;

switch(F.type)
 case 'gen'
  outsig=F.gs*insig;
  
 case 'dgt'
  outsig=idgt(framecoef2native(F,insig),F.gs,F.a);  
 case 'dgtreal'
  outsig=idgtreal(framecoef2native(F,insig),F.gs,F.a,F.M);  
 case 'dwilt'
  outsig=idwilt(framecoef2native(F,insig),F.gs);  
 case 'wmdct'
  outsig=iwmdct(framecoef2native(F,insig),F.gs);  
 
 case {'filterbank','ufilterbank'}
  outsig=ifilterbank(framecoef2native(F,insig),F.gs,F.a);   
 case {'filterbankreal','ufilterbankreal'}
  outsig=2*real(ifilterbank(framecoef2native(F,insig),F.gs,F.a));

 case {'nsdgt','unsdgt'}
  outsig=insdgt(framecoef2native(F,insig),F.gs,F.a);   
 case {'nsdgtreal','unsdgtreal'}
  outsig=insdgtreal(framecoef2native(F,insig),F.gs,F.a);
 
 case {'dcti','dctiv','dsti','dstiv'}
  outsig=feval(F.type,insig);
  
 case 'dctii'
  outsig=dctiii(insig);  
 case 'dctiii'
  outsig=dctii(insig);
 case 'dstii'
  outsig=dstiii(insig);
 case 'dstiii'
  outsig=dstii(insig);
 case 'dft'
  outsig=idft(insig);
 case 'fft'
  outsig=ifft(insig);
 case 'fftreal'
  outsig=ifftreal(insig,F.L);  
end;

  

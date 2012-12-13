function outsig=frsyn(F,insig);
%FRSYN  Frame synthesis operator
%   Usage: f=frsyn(F,c);
%
%   `f=frsyn(F,c)` constructs a signal *f* from the frame coefficients *c*
%   using the frame *F*. The frame object *F* must have been created using
%   |frame|_.
%
%   See also: frame, frana, plotframe
  
if nargin<2
  error('%s: Too few input parameters.',upper(mfilename));
end;

if ~isstruct(F)
  error('%s: First agument must be a frame definition structure.',upper(mfilename));
end;

switch(F.type)
  case 'identity'
    outsig=insig;    
  case 'gen'
    outsig=F.g*insig;
    
  case 'dgt'
    outsig=idgt(framecoef2native(F,insig),F.g,F.a,F.vars{:});  
  case 'dgtreal'
    outsig=idgtreal(framecoef2native(F,insig),F.g,F.a,F.M,F.vars{:});  
  case 'dwilt'
    outsig=idwilt(framecoef2native(F,insig),F.g);  
  case 'wmdct'
    outsig=iwmdct(framecoef2native(F,insig),F.g);  
    
  case {'filterbank','ufilterbank'}
    outsig=ifilterbank(framecoef2native(F,insig),F.g,F.a);   
  case {'filterbankreal','ufilterbankreal'}
    outsig=2*real(ifilterbank(framecoef2native(F,insig),F.g,F.a));

  case {'nsdgt','unsdgt'}
    outsig=insdgt(framecoef2native(F,insig),F.g,F.a);   
  case {'nsdgtreal','unsdgtreal'}
    outsig=insdgtreal(framecoef2native(F,insig),F.g,F.a);
    
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
  case 'fusion'
    W=size(insig,2);
    L=size(insig,1)/framered(F);

    outsig=zeros(L,W);

    idx=0;    
    for ii=1:F.Nframes
        coeflen=L*framered(F.frames{ii});
        outsig=outsig+frsyn(F.frames{ii},insig(idx+1:idx+coeflen,:))*F.w(ii);
        idx=idx+coeflen;
    end;
end;

  

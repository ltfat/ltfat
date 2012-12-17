function Fd=framedual(F);
%FRAMEDUAL  Construct the canonical dual frame
%   Usage: F=framedual(F);
%          F=framedual(F,L);
%
%   `Fd=frame(F)` returns the canonical dual frame of *F*.
%
%   See also: frame, framepair, frametight
  
if nargin<1
  error('%s: Too few input parameters.',upper(mfilename));
end;

% Default operation, work for a lot of frames
Fd=F;

% Handle the windowed transforms
switch(F.type)
  case {'dgt','dgtreal','dwilt','wmdct','filterbank','ufilterbank',...
        'nsdgt','unsdgt','nsdgtreal','unsdgtreal'}
    
    Fd.g={'dual',F.g};
    
  case {'filterbankreal','ufilterbankreal'}
    Fd.g={'realdual',F.g};  
    
  case 'gen'
    Fd.g=pinv(F.g)';
        
  case 'fusion'
    Fd.w=1./(F.Nframes*F.w);
    for ii=1:F.Nframes
        Fd.frames{ii}=framedual(F.frames{ii});
    end;
end;

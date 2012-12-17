function Ft=frametight(F);
%FRAMETIGHT  Construct the canonical tight frame
%   Usage: F=frametight(F);
%          F=frametight(F,L);
%
%   `Ft=frame(F)` returns the canonical tight frame of *F*.
%
%   See also: frame, framepair, framedual
  
if nargin<1
  error('%s: Too few input parameters.',upper(mfilename));
end;

% Default operation, work for a lot of frames
Ft=F;

% Handle the windowed transforms
switch(F.type)
  case {'dgt','dgtreal','dwilt','wmdct','filterbank','ufilterbank',...
        'nsdgt','unsdgt','nsdgtreal','unsdgtreal'}
    
    Ft.g={'tight',F.g};
    
  case {'filterbankreal','ufilterbankreal'}
    Ft.g={'realtight',F.g};  
    
  case 'gen'
    [U,sv,V] = svd(F.g,'econ');    
    Ft.g=U*V'; 
        
  case 'fusion'
    Ft.w=1./F.w;
    for ii=1:F.Nframes
        Ft.frames{ii}=frametight(F.frames{ii});
    end;
end;

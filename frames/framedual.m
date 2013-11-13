function Fd=framedual(F);
%FRAMEDUAL  Construct the canonical dual frame
%   Usage: F=framedual(F);
%          F=framedual(F,L);
%
%   `Fd=framedual(F)` returns the canonical dual frame of *F*.
%
%   The canonical dual frame can be used to get perfect reconstruction as in
%   the following example:::
%
%     % Create a frame and its canonical dual
%     F=frame('dgt','hamming',32,64);
%     Fd=framedual(F);
%
%     % Compute the frame coefficients and test for perfect
%     % reconstruction
%     f=gspi;
%     c=frana(F,f);
%     r=frsyn(Fd,c);
%     norm(r(1:length(f))-f)
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
        'nsdgt','unsdgt','nsdgtreal','unsdgtreal',...
        'fwt','ufwt','wfbt','uwfbt','wpfbt','uwpfbt'}
    
    Fd=frame(F.type,{'dual',F.g},F.origargs{2:end});
    
  case {'filterbankreal','ufilterbankreal'}
    Fd=frame(F.type,{'realdual',F.g},F.origargs{2:end});
    
  case 'gen'
    Fd=frame('gen',pinv(F.g)');
    
  case 'tensor'
    for ii=1:F.Nframes
        dual_frames{ii}=framedual(F.frames{ii});        
    end;
    Fd=frame('tensor',dual_frames{:});
        
  case 'fusion'
    dual_w=1./(F.Nframes*F.w);
    for ii=1:F.Nframes
        dual_frames{ii}=framedual(F.frames{ii});        
    end;
    Fd=frame('fusion',dual_w,dual_frames{:});
    
      
end;

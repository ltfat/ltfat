function Ft=frametight(F);
%FRAMETIGHT  Construct the canonical tight frame
%   Usage: F=frametight(F);
%          F=frametight(F,L);
%
%   `Ft=frametight(F)` returns the canonical tight frame of *F*.
%
%   The canonical tight frame can be used to get perfect reconstruction if
%   it is used for both analysis and synthesis. This is demonstrated in the
%   following example:::
%
%     % Create a frame and its canonical tight
%     F=frame('dgt','hamming',32,64);
%     Ft=frametight(F);
%
%     % Compute the frame coefficients and test for perfect
%     % reconstruction
%     f=gspi;
%     c=frana(Ft,f);
%     r=frsyn(Ft,c);
%     norm(r(1:length(f))-f)
%
%   See also: frame, framepair, framedual
  
if nargin<1
  error('%s: Too few input parameters.',upper(mfilename));
end;

% Default operation, works for a lot of frames
Ft=F;

% Handle the windowed transforms
switch(F.type)
  case {'dgt','dgtreal','dwilt','wmdct','filterbank','ufilterbank',...
        'nsdgt','unsdgt','nsdgtreal','unsdgtreal'}
    
    Ft=frame(F.type,{'tight',F.g},F.origargs{2:end});
    
  case {'filterbankreal','ufilterbankreal'}
    Ft=frame(F.type,{'realtight',F.g},F.origargs{2:end});
    
  case 'gen'
    [U,sv,V] = svd(F.g,'econ');    
    Ft=frame('gen',U*V');

  case 'tensor'
    for ii=1:F.Nframes
        tight_frames{ii}=frametight(F.frames{ii});
    end;
    F=frame('tensor',tight_frames{:});

  case 'fusion'
    tight_w=1./F.w;
    for ii=1:F.Nframes
        tight_frames{ii}=frametight(F.frames{ii});
    end;
    Ft=frame('fusion',tight_w,tight_frames{:});

end;

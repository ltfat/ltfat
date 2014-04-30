function Ft=frametight(F);
%FRAMETIGHT  Construct the canonical tight frame
%   Usage: F=frametight(F);
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
  
complainif_notenoughargs(nargin,1,'FRAMETIGHT');
complainif_notvalidframeobj(F,'FRAMETIGHT');

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
    
  case 'ufwt'
    % The canonical tight made from ufwt might not keep the iterated
    % filterbank structure
    [g,a] = wfbt2filterbank({F.g,F.J,'dwt'});
    g = comp_filterbankscale(g,a,F.flags.scaling);
    
    Ft = frametight(frame('ufilterbank',g,1,numel(g)));
    warning(sprintf(['%s: The canonical tight system does not preserve ',...
                     'the iterated filterbank structure.'],...
                     upper(mfilename)));
  case 'uwfbt'
    % The canonical tight made from ufwt might not keep the iterated
    % filterbank structure
    [g,a] = wfbt2filterbank(F.g,F.J);
    g = comp_filterbankscale(g,a,F.flags.scaling);
    
    Ft = frametight(frame('ufilterbank',g,1,numel(g)));
    warning(sprintf(['%s: The canonical tight system does not preserve ',...
                     'the iterated filterbank structure.'],...
                     upper(mfilename)));
  case 'uwpfbt'               
    % The canonical tight made from ufwt might not keep the iterated
    % filterbank structure
    [g, a] = wpfbt2filterbank(F.g,F.flags.interscaling);
    g = comp_filterbankscale(g,a,F.flags.scaling);
    
    Ft = frametight(frame('ufilterbank',g,1,numel(g)));
    warning(sprintf(['%s: The canonical tight system does not preserve ',...
                     'the iterated filterbank structure.'],...
                     upper(mfilename)));
   case 'wpfbt'
       error('TO DO:');
end;



% Treat the fixed length frames
if isfield(F,'fixedlength') && F.fixedlength
   Ft = frameaccel(Ft,F.L);
   Ft.fixedlength = 1;
end


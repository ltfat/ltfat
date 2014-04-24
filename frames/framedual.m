function Fd=framedual(F);
%FRAMEDUAL  Construct the canonical dual frame
%   Usage: F=framedual(F);
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
        'fwt','wfbt'}
    
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
    
  case 'ufwt'
    % The canonical dual of ufwt might not keep the iterated
    % filterbank structure
    [g, a] = wfbt2filterbank({F.g,F.J,'dwt'});
    g = comp_filterbankscale(g,a,F.flags.scaling);
    
    Fd = framedual(frame('ufilterbank',g,1,numel(g)));
    warning(sprintf(['%s: The canonical dual system does not preserve ',...
                     'the iterated filterbank structure.'],...
                     upper(mfilename)));
    
  case 'uwfbt'
    % The canonical dual of uwfbt might not keep the iterated
    % filterbank structure
    [g, a] = wfbt2filterbank(F.g);
    g = comp_filterbankscale(g,a,F.flags.scaling);
    
    Fd = framedual(frame('ufilterbank',g,1,numel(g)));
    warning(sprintf(['%s: The canonical dual system does not preserve ',...
                     'the iterated filterbank structure.'],...
                     upper(mfilename))); 
  case 'wpfbt'
    % The canonical dual of wpfbt might not keep the iterated
    % filterbank structure
    [g, a] = wpfbt2filterbank(F.g,F.flags.interscaling);      
    [gu,au,p] = nonu2ufilterbank(g,a);
    Fd = framedual(frame('ufilterbank',gu,au,numel(gu)));
    
    error('TO DO: ')
  case 'uwpfbt'
    % The canonical dual of uwfbt might not keep the iterated
    % filterbank structure
    [g, a] = wpfbt2filterbank(F.g,F.flags.interscaling);
    g = comp_filterbankscale(g,a,F.flags.scaling);
    
    Fd = framedual(frame('ufilterbank',g,1,numel(g)));
    warning(sprintf(['%s: The canonical dual system of frame type %s ',...
                     'does not preserve iterated filterbank structure.'],...
                     upper(mfilename),F.type));  
      
end;

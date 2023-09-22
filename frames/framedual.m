function Fd=framedual(F)
%FRAMEDUAL  Construct the canonical dual frame of a frame F
%   Usage: Fd=framedual(F);
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
  
complainif_notenoughargs(nargin,1,'FRAMEDUAL');
complainif_notvalidframeobj(F,'FRAMEDUAL');

% Default operation, work for a lot of frames
Fd=F;

% Handle the windowed transforms
switch(F.type)
  case {'dgt','dgtreal','dwilt','wmdct','filterbank','ufilterbank',...
        'nsdgt','unsdgt','nsdgtreal','unsdgtreal'}
    
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
    
    Fd = framedual(frame('filterbank',g,ones(numel(g),1),numel(g)));
%     warning(sprintf(['%s: The canonical dual system does not preserve ',...
%                      'the iterated filterbank structure.'],...
%                      upper(mfilename)));
    
  case 'uwfbt'
    % The canonical dual of uwfbt might not keep the iterated
    % filterbank structure
    [g, a] = wfbt2filterbank(F.g);
    g = comp_filterbankscale(g,a,F.flags.scaling);
    
    Fd = framedual(frame('filterbank',g,ones(numel(g),1),numel(g)));
%     warning(sprintf(['%s: The canonical dual system does not preserve ',...
%                      'the iterated filterbank structure.'],...
%                      upper(mfilename))); 
                 
  case 'uwpfbt'
    % The canonical dual of uwfbt might not keep the iterated
    % filterbank structure
    [g, a] = wpfbt2filterbank(F.g, F.flags.interscaling);
    g = comp_filterbankscale(g, a, F.flags.scaling);
    
    Fd = framedual(frame('filterbank',g,ones(numel(g),1),numel(g)));
%     warning(sprintf(['%s: The canonical dual system of frame type %s ',...
%                      'does not preserve iterated filterbank structure.'],...
%                      upper(mfilename),F.type));  
                 
  case 'fwt'
    is_basis = abs(sum(1./F.g.a)-1)<1e-6;
    is_tight = F.info.istight;
    
    % If the frame is a basis, there is only one dual frame.
    % If the basic filterbank is a parseval tight frame, the overal repr.
    % is also a tight frame.
    if is_basis || is_tight
        Fd = frame('fwt',{'dual',F.g},F.J);
    else
        error(['%s: Cannot create the canonical dual frame with the ',...
               'same structure. Consider casting the system to an ',...
               'uniform filterbank or using franaiter/frsyniter.'],...
               upper(mfilename)); 
    end

  case 'wfbt'
    is_basis = all(cellfun(@(nEl) abs(sum(1./nEl.a)-1)<1e-6,F.g.nodes));
    is_tight = F.info.istight;               
    
    if is_basis || is_tight
        Fd = frame('wfbt',{'dual',F.g});
    else
        error(['%s: Cannot create the canonical dual frame with the ',...
               'same structure. Consider casting the system to an ',...
               'uniform filterbank or using franaiter/frsyniter.'],...
               upper(mfilename)); 
    end
    
               
  case 'wpfbt'
     % WPFBT is too wierd.
     error(['%s: Canonical dual frame of wpfbt might not keep the ',...
           'same structure. Consider using franaiter/frsyniter.'],upper(mfilename));
    
end;

% Treat the fixed length frames
if isfield(F,'fixedlength') && F.fixedlength && isfield(F,'L')
   Fd = frameaccel(Fd,F.L);
   Fd.fixedlength = 1;
end




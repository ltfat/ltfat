function g = comp_filterbankscale(g,a,scaling)


if strcmp(scaling,'scale')
    g = cellfun(@(gEl,aEl) setfield(gEl,'h',gEl.h./aEl),g(:),...
                      num2cell(a(:)),...
                      'UniformOutput',0);
elseif strcmp(scaling,'noscale')
    % Do nothing
elseif strcmp(scaling,'sqrt')
    g = cellfun(@(gEl,aEl) setfield(gEl,'h',gEl.h./sqrt(aEl)),g(:),...
                      num2cell(a(:)),...
                      'UniformOutput',0); 
else
    error('%s: Unrecognized scaling flag.',upper(mfilename) );
end  


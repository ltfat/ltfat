function wt = comp_wpfbtscale(wt,interscaling)
%COMP_WPFBTSCALE Scale filters in the filterbank tree


% Here we handle scaling of intermediate outputs in the tree
if ~strcmpi(interscaling,'intnoscale')
    if strcmp('intscale',interscaling)
        interscalingfac = 1/2;
    elseif strcmp('intsqrt',interscaling)
        interscalingfac = 1/sqrt(2);
    end
    
    wtPath = nodeBForder(0,wt);
    rangeLoc = nodesLocOutRange(wtPath,wt);
    wtNodes = wt.nodes(wtPath);

    for ii=1:numel(wtPath)
          range = 1:numel(wtNodes{ii}.h);
          % Remove the outputs which are terminal
          range(rangeLoc{ii}) = [];
          wtNodes{ii}.h(range) = ...
             cellfun(@(hEl) setfield(hEl,'h',hEl.h*interscalingfac),...
                     wtNodes{ii}.h(range),...
                     'UniformOutput',0); 
                 
          wtNodes{ii}.g(range) = ...
             cellfun(@(hEl) setfield(hEl,'h',hEl.h*interscalingfac),...
                     wtNodes{ii}.g(range),...
                     'UniformOutput',0);
    end
    % Write the scaled ones back
    wt.nodes(wtPath) = wtNodes;
end


function c = comp_dtwfb(f,nodes,dualnodes,rangeLoc,rangeOut,ext,do_complex)

% First tree
c1 = comp_wfbt(f,nodes,rangeLoc,rangeOut,ext);
% Second tree
c2 = comp_wfbt(f,dualnodes,rangeLoc,rangeOut,ext);

% Combine outputs of trees
c = cellfun(@(crEl,ciEl) (crEl+1i*ciEl)/2,c1,c2,'UniformOutput',0);

if do_complex
   % Non-real specific
   cneg = cellfun(@(crEl,ciEl) (crEl-1i*ciEl)/2,c1,c2,...
                  'UniformOutput',0);

   c = [c;cneg(end:-1:1)];
end

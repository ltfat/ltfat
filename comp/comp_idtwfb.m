function f = comp_idtwfb(c,nodes,dualnodes,Lc,rangeLoc,rangeOut,ext,do_complex)


if do_complex
   % Split the coefficients
   c1 = cellfun(@(cEl1,cEl2) (cEl1 + cEl2)/2, c(1:end/2),c(end:-1:end/2+1),...
             'UniformOutput',0);
   c2 = cellfun(@(cEl1,cEl2) (-1i*cEl1 + 1i*cEl2)/2, c(1:end/2),c(end:-1:end/2+1),...
             'UniformOutput',0);
else
   c1 = cellfun(@real,c,'UniformOutput',0);
   c2 = cellfun(@imag,c,'UniformOutput',0);
end


f1 = comp_iwfbt(c1,nodes,Lc,rangeLoc,rangeOut,ext);
f = f1 + comp_iwfbt(c2,dualnodes,Lc,rangeLoc,rangeOut,ext); 

if ~do_complex
   f = real(f);
end
    

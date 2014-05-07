function L=comp_framelength_fusion(F,Ls);

% This is highly tricky: Get the minimal transform length for each
% subframe, and set the length as the lcm of that.
Lsmallest=1;
for ii=1:F.Nframes
    Lsmallest=lcm(Lsmallest,framelength(F.frames{ii},1));
end;
L=ceil(Ls/Lsmallest)*Lsmallest;

% Verify that we did not screw up the assumptions.
for ii=1:F.Nframes
    if L~=framelength(F.frames{ii},L)
        error(['%s: Cannot determine a frame length. Frame no. %i does ' ...
               'not support a length of L=%i.'],upper(mfilename),ii,L);
    end;
end;


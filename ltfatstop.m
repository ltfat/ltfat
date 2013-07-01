function ltfatstop()
%LTFATSTOP   Stops the LTFAT toolbox
%   Usage:  ltfatstop;
%
%   `ltfatstop` removes all LTFAT subdirectories from the path.
%
%   See also:  ltfatstart

%   AUTHOR : Peter L. Søndergaard.  

fullpath=strsplit(path,pathsep);

bp=ltfatbasepath;
% Remove the file separator at the end
bp=bp(1:end-1);
bplen=numel(bp);

for line=fullpath
    % Strip the cell array container away
    thispath=line{1};
    if numel(thispath)>=bplen && strcmp(thispath(1:bplen),bp)
        rmpath(thispath);
        disp(thispath)
    end;    
end;
    
  

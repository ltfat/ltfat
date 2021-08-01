function result=convert2Octave(result)

toOctave = 1;
if toOctave
  %for ii = 1:numel(result)
  for ii=1:1  
    %check if output is a cell array (like e.g. in filterbanks
    if iscell(result{ii})
      for jj=1:numel(result{ii})
        curres = result{ii}{jj};
        if isstruct(curres)
          curres = rewrite_struct2Octave(curres);
        end
        result{ii}{jj}=curres;
      end
    elseif isstruct(result{ii})
       curres = rewrite_struct2Octave(result{ii});
       result{ii}=curres;
    else
       curres = result{ii};
       if ischar(curres) && strcmp(curres(1), '@')
         result{ii}=str2func(curres);
       end
    end
  end 
end

%evalin('base', 'result')
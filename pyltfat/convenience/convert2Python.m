function result=convert2Python(result)

toOctave = 1;
if any(is_function_handle(result)) || any(isstruct(result)) || any(iscell(result))
  toOctave = 0;
  for ii = 1:numel(result)
    %check if output is a cell array (like e.g. in filterbanks
    if iscell(result{ii})
      %typeinfo(result{ii})
      for jj=1:numel(result{ii})
        curres = result{ii}{jj};
         %   typeinfo(result{ii}{jj})
        %if strcmp(typeinfo(result{ii}{jj}), 'sq_string')
        %  curres=char(curres)
        %end

        if isstruct(curres)
          [curres, toOctave] = rewrite_struct2Python(curres);
          %disp(typeinfo(curres.H))
        end
        result{ii}{jj}=curres;
      end
    elseif isstruct(result{ii})
       [curres, toOctave] = rewrite_struct2Python(result{ii});
       result{ii}=curres;     
    elseif is_function_handle(result{ii})
       result{ii} = char(result{ii});
       toOctave = 1;
    end

  end 
end
function _pyeval(input_file, output_file)
% _PYEVAL: Load a request from an input file, execute the request, and save
%         the response to the output file.
%
%   This allows you to run any Octave code. req should be a struct with the
%   following fields:
%       dname: The name of a directory to add to the runtime path before attempting to run the code.
%       func_name: The name of a function to invoke.
%       func_args: An array of arguments to send to the function.
%       nout: An int specifying how many output arguments are expected.
%       ref_indices: The indices of in the func_args that should
%         be replaced by the value represented by their name.
%       store_as: Optional name to store the return value in the base
%         workspace, instead of returing a value.
%
%   Should save a file containing the result object.
%
% Based on Max Jaderberg's web_feval
%
% This file is an extension of oct2py (Copyright (c) oct2py developers, distributed under the MIT License),
% to simplify calling the LTFAT from Python.
% Distributed under the terms of the GPL v3 License.

sentinel = { '__no_value__' };
result = { sentinel };
err = '';

try
    % Store the simple response in case we don't make it through the script.
    save('-v6', '-mat-binary', output_file, 'result', 'err');

    req = load(input_file);

    % Add function path to current path.
    if req.dname
        addpath(req.dname);
    end

    % Replace the names at the specified indices with their values.
    for idx=1:length(req.ref_indices)
      ref_index = req.ref_indices(idx);
      var_name = req.func_args{ref_index};
      req.func_args{ref_index} = evalin('base', var_name);
    end

    assignin('base', 'ans', sentinel);

    % Use the `ans` response if no output arguments are expected.
    if req.nout == 0

        if length(req.func_args)
          feval(req.func_name, req.func_args{:});
        else
          feval(req.func_name)
        end

        result = get_ans(sentinel);

    elseif length(req.func_args)
      try
        [result{1:req.nout}] = feval(req.func_name, req.func_args{:});
      catch ME
        if (strcmp(ME.message, 'element number 1 undefined in return list') != 1 ||
            length(ME.stack) != 1)
          rethrow(ME);
        else
          result = get_ans(sentinel);
        end

      end

    else
      try
        [result{1:req.nout}] = feval(req.func_name);
      catch ME
          rethrow(ME);
        if (strcmp(ME.message, 'element number 1 undefined in return list') != 1 ||
            length(ME.stack) != 1)
          rethrow(ME);
        end
      end
    end

    if req.store_as
      assignin('base', req.store_as, result{1});
      result = { sentinel };
    end

    if ((strcmp(get(0, 'defaultfigurevisible'), 'on') == 1) &&
        length(get(0, 'children')))
      drawnow('expose');
    end

catch ME
    err = ME;
end

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

%check if we have a function handle converted to a string
if 0 && toOctave
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
       curres = rewrite_struct2Octave(curres);
       result{ii}=curres;
    else
       curres = result{ii};
       if ischar(curres) && strcmp(curres(1), '@')
         result{ii}=str2func(curres);
       end
    end
  end 
end
%evalin('base', 'result{1}')
%result = [];
%check if output contains a function handle
%for ii = 1:length(result)
  
%    if isstruct(result{ii})
%        curres = result{ii};
%        fn = fieldnames(curres);
%        for k=1:numel(fn)
%            if( is_function_handle(curres.(fn{k})) )
%                curres.(fn{k}) = char(curres.(fn{k}));
%            end
%        end
%        result{ii} = curres;
%    else    
%        if is_function_handle(result{ii})
%            result{ii} = char(result{ii});
%        end
%    end    
%end
% Save the output to a file.
try
  save('-v6', '-mat-binary', output_file, 'result', 'err');
catch ME
  result = { sentinel };
  err = ME;
  save('-v6', '-mat-binary', output_file, 'result', 'err');
end

end  % function


function result = get_ans(sentinel)
    try
      [result{1}] = evalin('base', 'ans');
    catch
      result = { sentinel };
    end
end


function [out, ret] = rewrite_struct2Python(in)
    curres = in;
    ret = 0;
    fn = fieldnames(curres);
    for k=1:numel(fn)
        if( is_function_handle(curres.(fn{k})) )
            curres.(fn{k}) = func2str(curres.(fn{k}));
            ret = 1;%the function did something
        end
    end
    out = curres;
end

function out = rewrite_struct2Octave(in)
    curres = in;
    fn = fieldnames(curres);
    for k=1:numel(fn)
        if isfield(curres, fn(k)) && ischar(curres.(fn{k})) && strcmp(curres.(fn{k})(1), '@')
            curres.(fn{k}) = str2func(char(curres.(fn{k})));
        end
    end
    out = curres;
end
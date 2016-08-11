function g = filterbankscale(g,varargin)
%FILTERBANKSCALE Scale filters in filterbank
%   Usage:  g=filterbankscale(g,scal)
%           g=filterbankscale(g,'flag')
%           g=filterbankscale(g,L,'flag')
%
%   `g=filterbankscale(g,s)` scales each filter in *g* by multiplying it
%   with *scal*. *scal* can be either scalar or a vector of the same length
%   as *g*.
%   The function only works with filterbanks already instantiated 
%   (returned from a function with a `filters` suffix or run trough 
%   |filterbankwin|,|firfilter|,|blfilter|) such that the elements
%   of *g* must be either structs with .h or .H fields or be plain 
%   numeric vectors.
%
%   `g=filterbankscale(g,'flag')` instead normalizes each filter to have 
%   unit norm defined by 'flag'. It can be any of the flags recognized by 
%   |normalize|. The  normalization is done in the time domain by default. 
%   The normalization can be done in frequency by passing extra flag 'freq'.
%
%   `g=filterbankscale(g,L,'flag')` works as before, but some filters require 
%   knowing *L* to be instantialized to obtain their norm. The normalization
%   will be valid for the lengh *L* only. 
%
%   In any case, the returned filters will be in exactly the same format as 
%   the input filters.
%

complainif_notenoughargs(nargin,2,'FILTERBANKSCALE');

definput.import={'normalize'};
definput.importdefaults={'norm_notset'};
definput.flags.normfreq = {'nofreq','freq'};
definput.keyvals.arg1 = [];
[flags,kv,arg1]=ltfatarghelper({'arg1'},definput,varargin);

% Try running filterbankwin without L. This should fail
% for any strange filter definitions like 'gauss','hann',{'dual',...}
try
   g2 = filterbankwin(g,1,'normal');
catch
   err = lasterror;
   if strcmp(err.identifier,'L:undefined')
       % If it blotched because of the undefined L, explain that.
       % This should capture only formats like {'dual',...} and {'gauss'}
       error(['%s: Function cannot handle g in such format. ',...
              'Consider pre-formatting the filterbank by ',...
              'calling g = FILTERBANKWIN(g,a) or ',...
              'g = FILTERBANKWIN(g,a,L) first.'],upper(mfilename));
   else
       % Otherwise just rethrow the error
       error(err.message);
   end
end


% At this point, elements of g can only be:
% struct with numeric field .h,
% struct with numeric field .H
% struct with function handle in .H
% numeric vectors

if flags.do_norm_notset
   % No flag from normalize was set
   s = scalardistribute(arg1,ones(size(g)));

   for ii=1:numel(g)
       if isstruct(g{ii})
           % Only work with .h or .H, any other struct field is not
           % relevant
           if isfield(g{ii},'h')
               if ~isnumeric(g{ii}.h)
                   error('%s: g{ii}.h must be numeric',upper(mfilename));
               end
               g{ii}.h = s(ii)*g{ii}.h;
           elseif isfield(g{ii},'H')
               if isa(g{ii}.H,'function_handle')
                   g{ii}.H = @(L) s(ii)*g{ii}.H(L);
               elseif isnumeric(g{ii}.H)
                   g{ii}.H = s(ii)*g{ii}.H;
               else
                   error(['%s: g{ii}.H must be either numeric or a ',...
                          ' function handle'],upper(mfilename));
               end
           else
              error('%s: SENTINEL. Unrecognized filter struct format',...
                    upper(mfilename));
           end
       elseif isnumeric(g{ii})
           % This is easy
           g{ii} = s(ii)*g{ii};
       else
           error('%s: SENTINEL. Unrecognized filter format',...
                 upper(mfilename));
       end
   end
else
   % Normalize flag was set
   for ii=1:numel(g)
      L = arg1; % can be still empty
      % Run again with L specified
      [g2,~,info] = filterbankwin(g,1,L,'normal');

      if ~isempty(L) && L < max(info.gl)
          error('%s: One of the windows is longer than the transform length.',upper(mfilename));
      end;

      if isstruct(g{ii})
          if isfield(g{ii},'h')
              if ~isnumeric(g{ii}.h)
                  error('%s: g{ii}.h must be numeric',upper(mfilename));
              end
              % Normalize either in time or in the frequency domain

              if flags.do_freq
                  complain_L(L);
                  % Get frequency response and it's norm
                  H = comp_transferfunction(g2{ii},L);
                  [~,fn] = normalize(H,flags.norm);
                  g{ii}.h = g{ii}.h/fn;
              else
                  if isfield(g{ii},'fc') && g{ii}.fc~=0
                      complain_L(L); % L is required to do a proper modulation
                  else
                      L = numel(g2{ii}.h);
                  end
                  % Get impulse response with all the fields applied
                  tmpg = comp_filterbank_pre(g2(ii),1,L,inf);
                  [~,fn] = normalize(tmpg{1}.h,flags.norm);
                  g{ii}.h = g{ii}.h/fn; 
              end
          elseif isfield(g{ii},'H')
              if isa(g{ii}.H,'function_handle')
                  complain_L(L);
                  H = comp_transferfunction(g2{ii},L);
                  if flags.do_freq
                     [~,fn] = normalize(H,flags.norm);
                     g{ii}.H = @(L) g{ii}.H(L)/fn;
                  else
                     [~,fn] = normalize(ifft(H),flags.norm);
                     g{ii}.H = @(L) g{ii}.H(L)/fn;
                  end
              elseif isnumeric(g{ii}.H)
                  if ~isfield(g{ii},'L')
                      error('%s: g.H is numeric, but .L field is missing',...
                            upper(mfilename));
                  end
                  if isempty(L)
                      L = g{ii}.L;
                  else
                      if L ~= g{ii}.L
                          error('%s: L and g.L are not equal',...
                                upper(mfilename));
                      end
                  end
                  
                  H = comp_transferfunction(g2{ii},L);
                  if flags.do_freq
                     [~,fn] = normalize(H,flags.norm);
                     g{ii}.H = g{ii}.H/fn;
                  else
                     [~,fn] = normalize(ifft(H),flags.norm);
                     g{ii}.H = g{ii}.H/fn;
                  end
              end
          else
             error('%s: SENTINEL. Unrecognized filter struct format',...
                   upper(mfilename));
          end
      elseif isnumeric(g{ii})
          % This one is not so easy
          if flags.do_freq
             complain_L(L);
             % We must use g2 here
             [~, fn] = normalize(fft(g2{ii}.h,L),flags.norm);
             g{ii} = g{ii}/fn;
          else
             g{ii} = normalize(g{ii},flags.norm); 
          end
      else
          error('%s: SENTINEL. Unrecognized filter format',...
                upper(mfilename));
      end
   end
end

function complain_L(L)

if isempty(L)
     error('%s: L must be specified',upper(mfilename));
end

complainif_notposint(L,'L',mfilename)

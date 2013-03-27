function [w] = fwtinit(wdef,varargin)
%FWTINIT  Fast Wavelet Transform Filterbank Structure Initialization
%   Usage:  w = fwtinit(wdef,...);
%           w = fwtinit(wdef,...,'a',a);
%
%   Input parameters:
%         wdef : Wavelet filters specification.
%
%   Output parameters:
%         w       : Structure defining the filters.
%
%   
%   `w=fwtinit(wdef,...)` produces a structure describing the analysis
%   and synthesis parts of a one level perfect reconstruction wavelet-type
%   filterbank. The function is a wrapper for calling all the functions 
%   starting with `wfilt_` defined in the LTFAT wavelets directory.
%   The possible formats of the `wdef` are the following:
%
%     1) One option is passing a cell array whose first element is the
%        name of the function defining the basic wavelet filters (`wfilt_`
%        prefix) and the other elements are the parameters of the
%        function. e.g. `{'db',10}` calls `wfilt_db(10)` internally.
%
%     2) Character string as concatenation of the name of the wavelet
%        filters defining function (as above) and the numeric parameters
%        delimited by ':' character, e.g. 'db10' has the same effect as above,
%        'spline4:4' calls `wfilt_spline(4,4)` internally.
%
%     3) The third possible format of $h$ is to pass cell array of one
%        dimensional numerical vectors directly defining the wavelet filter
%        impulse responses.  By default, outputs of the filters are
%        subsampled by a factor equal to the number of the filters. Pass
%        additional scalar or vector *a* to define custom subsampling
%        factors.
%
%     4) The fourth option is to pass a structure obtained from the
%        |fwtinit| function. The function check if the struct is valid
%        and update `w.filts` and `w.a` fields according to the other
%        parameters passed.
%
%   The output structure has the following fields:
%
%     `w.filts`
%        analysis or synthesis filterbank
%
%     `w.h`
%        analysis filter bank
%     
%     `w.g`
%        synthesis filter bank
%
%     `w.a`
%        implicit subsampling factors
% 
%   Remark: Function names starting with `wfilt_` cannot contain numbers! 
%
%   Choosing wavelet filters
%   ------------------------
%   
%   Determining which wavelet filters to use depends strongly on the type
%   of the analyzed signal. There are several properties which should be
%   taken into account.
%   
%   orthogonality/biorthogonality/frame
%   number of vanishing moments psi
%   symmetry/linear phase
%
%   support of psi
%   smoothnes and regularity of psi
%   
%   
%
%   See also: fwt, ifwt, wfilt_db
%
%   References: ma08wt


% chached last filterbank
persistent cachw;
% cached function parameters passed last
persistent cachwDesc;

% wavelet filters functions definition prefix
wprefix = 'wfilt_';
waveletsDir = 'wavelets';
numDelimiter = ':';


% output structure definition
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
w.filts = {};
w.h = {};
w.g = {};
w.a = [];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% return empty struct if no argument was passed
if nargin<1
  return;
end;

% process other parameters

definput.import = {'fwtcommon'};
% Contents of the arg_fwtcommon:
% definput.flags.ansy = {'ana','syn'};
% definput.keyvals.a = [];
[flags,kv]=ltfatarghelper({},definput,varargin);

% Was the function called before with the same parameters?
% if yes, return the chached one
if(isequal(cachwDesc,{wdef,kv.a}))
   w = cachw;
   w = updateTransDirect(flags.do_ana,w);
   return;
else
   cachwDesc = {wdef,kv.a};
end

if ischar(wdef)
   % Process wdef in format 2)
   try
    wname = parseNameValPair(wdef,wprefix);
   catch err
      % If failed, clean the cache.
      cachwDesc = [];
      cachw = [];
      error(err.message);
   end
elseif iscell(wdef)
    % Process wdef in format 3)
    wname = wdef;
    if ~ischar(wname{1})
       if iscell(wname{1})
          if(length(wname)==1)
              w = formatFilters(wname,flags.do_ana,[],w);
              w = updateTransDirect(flags.do_ana,w);
              cachw = w;
              return;
          else
             error('%s: Unrecognizer format of the filterbank definition.',upper(mfilename)); 
          end
       elseif isnumeric(wname{1})
            w = formatFilters(wname,flags.do_ana,[],w);
            w = updateTransDirect(flags.do_ana,w);
       else
          error('%s: Unrecognizer format of the filterbank definition.',upper(mfilename));
       end
       w.a = formata(max([length(w.h),length(w.g)]),kv.a,w.a);
       cachw = w;
       return;
    end
elseif isstruct(wdef)
    % Process wdef in format 4)
    if(isequal(fieldnames(wdef),fieldnames(w)))
        w = wdef;
        w = updateTransDirect(flags.do_ana,w);
        w.a = formata(length(w.h),kv.a,w.a);
        cachw = w;
        return;
    else
       error('%s: Passed structure has different fields.',upper(mfilename)); 
    end
else
    error('%s: First argument must be a string, cell or struct.',upper(mfilename));
end;

% wname now contains wdef in format 1)

% Search for m-file containing string wname
wfiltFile = dir(fullfile(ltfatbasepath,sprintf('%s/%s%s.m',waveletsDir,wprefix,lower(wname{1}))));
if(isempty(wfiltFile))
   error('%s: Unknown wavelet type: %s',upper(mfilename),name); 
else
   % if found, crop '.m' from the filename 
   tmpFile = wfiltFile.name(1:end-2); 
end

% Synthesis filters delay
d = [];
wfiltNargout = nargout(tmpFile);

if(nargin(tmpFile)~=numel(wname)-1)
   error('%s: Incorrect number of parameters to be passed to the wfilt_ func.',upper(mfilename));
end

if(wfiltNargout==3)
   [tmph, tmpg, w.a] = feval(tmpFile,wname{2:end});
elseif(wfiltNargout==4) 
   [tmph, tmpg, w.a, d] = feval(tmpFile,wname{2:end});
else
   error('%s: Function %s does not return 3 or 4 arguments.',upper(mfilename),upper(tmpFile));
end

w.a = formata(length(w.h),kv.a,w.a);

w = formatFilters(tmph,1,d,w);
w = formatFilters(tmpg,0,d,w);
w = updateTransDirect(flags.do_ana,w);

cachw = w;

function w = formatFilters(cellh,do_ana,d,w)
   noFilts = numel(cellh);
   if(isempty(d))
      d = findFiltDelays(cellh,do_ana,'half');
   end
   if(do_ana)
      w.h = cell(noFilts,1);
      for ff=1:noFilts
         w.h{ff} = wfiltstruct('FIR');
         w.h{ff}.h = cellh{ff};
         w.h{ff}.d = d(ff);
      end
   else
      w.g = cell(noFilts,1);
      for ff=1:noFilts
         w.g{ff} = wfiltstruct('FIR');
         w.g{ff}.h = cellh{ff};
         w.g{ff}.d = d(ff);
      end 
   end
end %END FORMATFILTERS

function w = updateTransDirect(do_ana,w)
   if(do_ana)
       w.filts = w.h;
   else
       w.filts = w.g;
   end
end %END UPDATETRANSDIRECT


function a = formata(filtsNo,a,origa)
   if(~isempty(a))
      if(length(a)==1)
         atmp = ones(filtsNo,1);
         a = atmp*a;
      end
   else
       if(~isempty(origa))
          a = origa;
       else
          atmp = ones(filtsNo,1);
          a = atmp*filtsNo;
       end
   end
end %END FORMATA

function wcell = parseNameValPair(wchar,wprefix)
%PARSENAMEVALPAIR
%Parses string in the following format wnameN1:N2... , where wname have to
%be name of the existing function with wfilt_ prefix. N1,N2,... are doubles
%delimited by character ':'.
%The output is cell array {wname,str2double(N1),str2double(N2),...}
%The wfilt_ function name cannot contain numbers

wcharNoNum = wchar(1:find(isstrprop(wchar,'digit')~=0,1)-1);

wfiltFiles = dir(fullfile(ltfatbasepath,sprintf('%s/%s*',waveletsDir,wprefix)));
wfiltNames = arrayfun(@(fEl) fEl.name(1+find(fEl.name=='_',1):find(fEl.name=='.',1,'last')-1),wfiltFiles,'UniformOutput',0);
wcharMatch = cellfun(@(nEl) strcmpi(wcharNoNum,nEl),wfiltNames);
wcharMatchIdx = find(wcharMatch~=0);
if(isempty(wcharMatchIdx))
   error('%s: Unknown wavelet filter definition string.',upper(mfilename));
end
if(numel(wcharMatchIdx)>1)
   error('%s: Ambiguous wavelet filter definition string. Probably bug somewhere.',upper(mfilename));
end

wcell{1} = wfiltNames{wcharMatchIdx};
numString = wchar(numel(wcell{1})+1:end);
if(isempty(numString))
   error('%s: No numeric parameter specified in %s.',upper(mfilename),wchar); 
end
wcharNum = textscan(numString,'%f','Delimiter',numDelimiter);
if(~isnumeric(wcharNum{1})||any(isnan(wcharNum{1})))
   error('%s: Incorrect numeric part of the wavelet filter definition string.',upper(mfilename));
end
wcell = [wcell, num2cell(wcharNum{1}).'];
end %END PARSENAMEVALPAIR

function d = findFiltDelays(cellh,do_ana,type)
   filtNo = numel(cellh);
   d = ones(filtNo,1);

   for ff=1:filtNo
       if(strcmp(type,'half'))
           if(do_ana)
               d(ff) = ceil((length(cellh{ff})+1)/2);
           else
               d(ff) = floor((length(cellh{ff})+1)/2);
           end
       elseif(strcmp(type,'zeroana'))
            if(do_ana)
               d(ff) = 1;
           else
               d(ff) = length(cellh{ff});
            end
       elseif(strcmp(type,'zerosyn'))
            if(do_ana)
               d(ff) = length(cellh{ff});
           else
               d(ff) = 1;
            end
       elseif(strcmp(type,'energycent'))
           tmphh =cellh{ff};
           tmphLen = length(tmphh);
           ecent = sum((1:tmphLen-1).*tmphh(2:end).^2)/sum(tmphh.^2);
           if(do_ana)
               d(ff) = round(ecent)+1;
               if(rem(abs(d(ff)-d(1)),2)~=0)
                  d(ff)=d(ff)+1;
               end
           else
               anad = round(ecent)+1;
               d(ff) = tmphLen-anad;
               if(rem(abs(d(ff)-d(1)),2)~=0)
                  d(ff)=d(ff)-1;
               end
           end

       else        
           error('fail');
       end
  end
end %END FINDFILTDELAYS

end %END FWTINIT







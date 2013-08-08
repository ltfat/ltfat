function [w] = fwtinit(wdef)
%FWTINIT  Wavelet Filterbank Structure Initialization
%   Usage:  w = fwtinit(wdef);
%
%   Input parameters:
%         wdef : Wavelet filters specification.
%
%   Output parameters:
%         w    : Structure defining the filterbank.
%   
%   `w=fwtinit(wdef)` produces a structure describing the analysis 
%   (field `.h`) synthesis (field `.g`) parts of a one level wavelet-type 
%   filterbank. The function is a wrapper for calling all the functions 
%   starting with `wfilt_` defined in the LTFAT wavelets directory.
%
%   The possible formats of the `wdef` are the following:
%
%   1) Cell array whose first element is the name of the function defining
%      the basic wavelet filters (`wfilt_` prefix) and the other elements
%      are the parameters of the function. e.g. `{'db',10}` calls 
%      `wfilt_db(10)` internally.
%
%   2) Character string as concatenation of the name of the wavelet
%      filters defining function (as above) and the numeric parameters
%      delimited by ':' character, e.g. `'db10'` has the same effect as above,
%      `'spline4:4'` calls `wfilt_spline(4,4)` internally.
%
%   3) Cell array of one dimensional numerical vectors directly defining
%      the wavelet filter impulse responses.  By default, outputs of the 
%      filters are subsampled by a factor equal to the number of the 
%      filters. Pass additional key-value pair 'a',a (still inside of the
%      cell array) to define the custom subsampling factors, e.g.: 
%      {h1,h2,'a',[2,2]}.
%
%   4) The fourth option is to pass again the structure obtained from the
%      |fwtinit| function. The structure is checked whether it has a valid
%      format.
%
%   5) Two element cell array. First element is the string `'dual'` and the
%      second one is in format 1), 2) or 4). This returns dual of whatever
%      is passed.
%
%   6) Two element cell array. First element is the string `'strict'` and the
%      second one is in format 1), 2), 4) or 5). This ensures the wavelet
%      filters not forming a tight frame has to be explicitly defined.
%
%   Using synthesis filters for analysis and vice versa makes a difference
%   in case of birtoghonal filters (e.g. `'spline4:4'`) and filters which 
%   constructs a general frame (e.g. `'symds2'`), in other cases, the analysis
%   filters are just time-reversed version of synthesis filters. To specify
%   which filters should be used, one can re-use the items 1) and 2) in the
%   following way:
%
%   1) Add `'ana'` or `'syn'` as the first element in the cell array e.g. 
%      `{'ana','spline',4,4}` or `{'syn','spline',4,4}`.
%
%   2) Add `'ana:'` or `'syn:'` to the beginning of the string e.g. 
%      `'ana:spline4:4'` or `'syn:spline4:4'`.
%
%   Please note that using e.g. `c=fwt(f,'ana:spline4:4',J)` and 
%   `fhat=ifwt(c,'ana:spline4:4',J,size(f,1))` will not give a perfect
%   reconstruction.
%
%   The output structure has the following fields:
%
%     w.h  Cell array of structures defining the original analysis filter bank.
%     
%     w.g  Cell array of structures defining the original synthesis filter bank.
%
%     w.a  Subsampling factors.
%
%     w.origArgs 
%          Original parameters in format 1).
%
%   `w.h`, `w.g` are cell-arrays of structures. Each structure represents 
%   one filter in the filterbank and it has the following fields:
%
%     .offset  Filter delay in samples.
%
%     .h       Actual filter impulse response values.
%
%   **Remark:** Function names with the `wfilt_` prefix cannot contain numbers
%   and cannot start with 'ana' or 'syn'! 
%
%   Choosing wavelet filters
%   ------------------------
%   
%   TO DO: Write something meaningfull.
%   Determining which wavelet filters to use depends strongly on the type
%   of the analyzed signal. There are several properties which should be
%   taken into account.
%   
%   orthogonality/biorthogonality/frame
%   number of vanishing moments of psi
%   symmetry/linear phase
%
%   support of psi
%   smoothnes/regularity of psi
%   
%
%   See also: fwt, ifwt, wfilt_db
%
%   References: ma08wt


% chached last filterbank
% persistent cachw;
% cached function parameters passed last
% persistent cachwDesc;



% wavelet filters functions definition prefix
wprefix = 'wfilt_';
waveletsDir = 'wavelets';


% output structure definition
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
w.origArgs = {};
w.h = {};
w.g = {};
w.a = [];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% return empty struct if no argument was passed
if nargin<1
  return;
end;

if isempty(wdef)
    error('%s: Input argument is empty.',upper(mfilename)); 
end


do_strict = 0;
do_dual = 0;

% Check 'strict'
if iscell(wdef) && ischar(wdef{1}) && strcmpi(wdef{1},'strict')
   do_strict = 1;
   wdef = wdef{2:end};
end
if iscell(wdef) && ischar(wdef{1}) && strcmpi(wdef{1},'dual')
   do_dual = 1;
   wdef = wdef{2:end};
end

% Disabled caching
% Was the function called before with the same parameters?
% if yes, return the chached one
%  if(isequal(cachwDesc,wdef))
%     complaina(kv.a,'structure');
%     w = cachw;
%     w.filts = decideFilters(w.forceAna,ansyBool,w);
%     return;
%  else
%     cachwDesc = wdef;
%  end

if isstruct(wdef)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Process wdef in format 4)%
    % Do checks and return quicky %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Check the fields
    if isequal(fieldnames(wdef),fieldnames(w))
        if ~do_dual && ~do_strict
           w = wdef;
           %cachw = w;
           return;
        else 
           if ~isempty(wdef.origArgs)
              wdef = wdef.origArgs;
           else
              error('%s: The structure was not buit using compatible formats.',upper(mfilename));
           end
        end
    else
       error('%s: Passed structure has different fields.',upper(mfilename)); 
    end
end

if iscell(wdef)
    wname = wdef;
    if ~ischar(wname{1})
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%
       % Process wdef in format 3)%
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%

       if isnumeric(wname{1})
             complainDual(do_dual,'numeric cell array');
             equalsa = cellfun(@(wEl)strcmp(wEl,'a'),wname);
             apos = find(equalsa==1);
             if isempty(apos)
                apos = numel(wname)+1;
                w.a = ones(numel(wname),1)*numel(wname);
             else
                if apos==numel(wname)-1 && isnumeric(wname{apos+1}) && numel(wname{apos+1})==apos-1
                   w.a = wname{apos+1};
                else
                   error('%s: Key ''a'' have to be followed by a vector of length %i.',upper(mfilename),apos-1);
                end
             end
             w.h = formatFilters(wname(1:apos-1),[]);
             w.g = formatFilters(wname(1:apos-1),[]);
       else
          error('%s: Unrecognizer format of the filterbank definition.',upper(mfilename));
       end
       
       %cachw = w;
       return;
    end
elseif ischar(wdef)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Process wdef in format 2)%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
   try
       wname = parseNameValPair(wdef,wprefix);
       % Octave does not support the "catch err" stament, so use "lasterror"
       % instead    
       %catch err
   catch
      err=lasterror;
      % If failed, clean the cache.
      cachwDesc = [];
      cachw = [];
      error(err.message);
   end
else
    error('%s: First argument must be a string, cell or struct.',upper(mfilename));
end;

do_forceAna = [];
is_tight = 0;

% Check whether wavelet definition starts with ana or syn
if ischar(wname{1}) && numel(wname{1})==3
   if strcmpi(wname{1},'ana') || strcmpi(wname{1},'syn')
      % Set field only if ana or syn was explicitly specified.
      do_forceAna = strcmpi(wname{1},'ana');
      wname = wname(2:end);
   end
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% wname now contains wdef in format 1)%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Search for m-file containing string wname
wfiltFile = dir(fullfile(ltfatbasepath,sprintf('%s/%s%s.m',waveletsDir,wprefix,lower(wname{1}))));
if(isempty(wfiltFile))
   error('%s: Unknown wavelet type: %s',upper(mfilename),wname{1}); 
else
   % if found, crop '.m' from the filename 
   tmpFile = wfiltFile.name(1:end-2); 
end

% There is a bug in nargout in version 3.6 of Octave, but not in later
% stable versions
if isoctave
    octs=strsplit(version,'.');
    octN=str2num(octs{1})*1000+str2num(octs{2});
    if octN<3008
        try
            feval(tmpFile);
        catch
        end;
    end;
end;
    
wfiltNargout = nargout(tmpFile);

if(nargin(tmpFile)~=numel(wname)-1)
   error('%s: Incorrect number of parameters to be passed to the %s func.',upper(mfilename),tmpFile);
end

info = [];
if(wfiltNargout==3)
   [tmph, tmpg, w.a] = feval(tmpFile,wname{2:end});
elseif(wfiltNargout==4) 
   [tmph, tmpg, w.a, info] = feval(tmpFile,wname{2:end});
else
   error('%s: Function %s does not return 3 or 4 arguments.',upper(mfilename),upper(tmpFile));
end


if ~isempty(info)&&isfield(info,'istight')
   is_tight = info.istight;
end

d = [];
if isfield(info,'d')
   d = info.d;
end

if numel(tmph)~=numel(w.a) || numel(tmpg)~=numel(w.a)
   error('%s: Variables returned by %s have different element counts.',upper(mfilename),upper(tmpFile));
end

if ~is_tight && do_strict && isempty(do_forceAna)
   error(['%s: %s filters does not form a tight frame. Choose either ''ana:%s'' or ''syn:%s'' '],upper(mfilename),tmpFile,wcell2str(wname),wcell2str(wname));
end
 
w.h = formatFilters(tmph,d);
w.g = formatFilters(tmpg,d);

w.origArgs = wname;

if ~isempty(do_forceAna)
   if do_dual
      do_forceAna = ~do_forceAna; 
   end
   if do_forceAna
      w.g = w.h;
      w.origArgs = {'ana',w.origArgs{:}};
   else
      w.h = w.g;
      w.origArgs = {'syn',w.origArgs{:}};
   end
end


end %END FWTINIT

function filts = formatFilters(cellf,d)
   noFilts = numel(cellf);
   filts = cell(noFilts,1);
   if(isempty(d))
      d = findFiltDelays(cellf,'half');
   end

   for ff=1:noFilts
      %filts{ff} = wfiltstruct('FIR');
      filts{ff}.h = cellf{ff}(:);
      filts{ff}.offset = -d(ff);
   end

end %END FORMATFILTERS

function wcell = parseNameValPair(wchar,wprefix)
%PARSENAMEVALPAIR
%Parses string in the following format wnameN1:N2... , where wname have to
%be name of the existing function with wfilt_ prefix. N1,N2,... are doubles
%delimited by character ':'.
%The output is cell array {wname,str2double(N1),str2double(N2),...}
%The wfilt_ function name cannot contain numbers

wcell = {}; 
numDelimiter = ':';

% Check whether the first 4 characters are 'ana:' or 'syn:'
if numel(wchar)>4
   if strcmpi(wchar(1:4),'ana:')
      wcell = [wcell,{'ana'}];
      wchar = wchar(5:end);
   elseif strcmpi(wchar(1:4),'syn:')
      wcell = [wcell,{'syn'}];
      wchar = wchar(5:end);
   end
end

% Take out all numbers from the string
wcharNoNum = wchar(1:find(isstrprop(wchar,'digit')~=0,1)-1);

% List all files satysfying the following: [ltfatbase]/wavelets/wfilt_*.m?
wfiltFiles = dir(fullfile(ltfatbasepath,sprintf('%s/%s*.m','wavelets',wprefix)));
% Get just the filanames without the wfilt_ prefix
wfiltNames = arrayfun(@(fEl) fEl.name(1+find(fEl.name=='_',1):find(fEl.name=='.',1,'last')-1),wfiltFiles,'UniformOutput',0);
% Compare the filenames with a given string
wcharMatch = cellfun(@(nEl) strcmpi(wcharNoNum,nEl),wfiltNames);
% Find index(es) of the matches.
wcharMatchIdx = find(wcharMatch~=0);
% Handle faulty results.
if(isempty(wcharMatchIdx))
   dirListStr = cell2mat(cellfun(@(wEl) sprintf('%s, ',wEl), wfiltNames(:)','UniformOutput',0));
   error('%s: Unknown wavelet filter definition string: %s.\nAccepted are:\n%s',upper(mfilename),wcharNoNum,dirListStr(1:end-2));
end
if(numel(wcharMatchIdx)>1)
   error('%s: Ambiguous wavelet filter definition string. Probably bug somewhere.',upper(mfilename));
end


match = wfiltNames{wcharMatchIdx};
wcell = [wcell,{match}];
% Extract the numerical parameters from the string (delimited by :)
numString = wchar(numel(match)+1:end);
if(isempty(numString))
   error('%s: No numeric parameter specified in %s.',upper(mfilename),wchar); 
end
% Parse the numbers.
wcharNum = textscan(numString,'%f','Delimiter',numDelimiter);
if(~isnumeric(wcharNum{1})||any(isnan(wcharNum{1})))
   error('%s: Incorrect numeric part of the wavelet filter definition string.',upper(mfilename));
end
wcell = [wcell, num2cell(wcharNum{1}).'];
end %END PARSENAMEVALPAIR

function d = findFiltDelays(cellh,type)
   filtNo = numel(cellh);
   d = ones(filtNo,1);

   for ff=1:filtNo
       if(strcmp(type,'half'))
               d(ff) = floor((length(cellh{ff})+1)/2);
%        elseif(strcmp(type,'energycent'))
%            tmphh =cellh{ff};
%            tmphLen = length(tmphh);
%            ecent = sum((1:tmphLen-1).*tmphh(2:end).^2)/sum(tmphh.^2);
%            if(do_ana)
%                d(ff) = round(ecent)+1;
%                if(rem(abs(d(ff)-d(1)),2)~=0)
%                   d(ff)=d(ff)+1;
%                end
%            else
%                anad = round(ecent)+1;
%                d(ff) = tmphLen-anad;
%                if(rem(abs(d(ff)-d(1)),2)~=0)
%                   d(ff)=d(ff)-1;
%                end
%            end

       else        
           error('TO DO: Unsupported type.');
       end
  end
end %END FINDFILTDELAYS

function complainDual(dual,whereStr)
if dual
  error('%s: ''dual'' option not allowed for the %s input.',upper(mfilename),whereStr); 
end
end % END COMPLAINA

function str = wcell2str(wcell)
strNums = cellfun(@(wEl) [num2str(wEl),':'],wcell(2:end),'UniformOutput',0);
strNums = cell2mat(strNums);
str = [wcell{1},strNums(1:end-1)];
end








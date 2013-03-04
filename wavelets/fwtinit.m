function [w] = fwtinit(wavname,varargin)
%FWTINIT  Fast Wavelet Transform Filterbank Structure Initialization
%   Usage:  w = fwtinit(wavname,'a',a);
%
%   Input parameters:
%         wavname : Wavelet filters generating function name (without the prefix).
%
%   Output parameters:
%         w       : Structure defining the filters.
%
%   `w=fwtinit(wavname,...)` produces a structure describing the analysis
%   and synthesis parts of a one level perfect reconstruction wavelet-type
%   filterbank.  The *wavname* can be a string defining a function name (without
%   the prefix) to be called or it can be a cell-array with the first element
%   beeing the function name and the remaining elements are passed further
%   to the `wfilt_` function.
%
%   The function is a wrapper for calling all the functions starting with
%   `wfilt_` defined in the LTFAT wavelets directory. The structure which
%   the function produces can be directly passed to all
%   functions instead of the cell-arrays with wavelet filters which all
%   `wfilt_` functions produces.
%
%   The structure has the following fields:
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
%
%   Choosing wavelet filters
%   ------------------------
%   
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

% was the function called before with the same parameters?
% if yes, return the chached one
if(isequal(cachwDesc,{wavname,kv.a}))
   w = cachw;
   w = updateTransDirect(flags.do_ana,w);
   return;
else
   cachwDesc = {wavname,kv.a};
end

if ischar(wavname)
    wname = {wavname};
elseif iscell(wavname)
    wname = wavname;
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
elseif isstruct(wavname)
    if(isequal(fieldnames(wavname),fieldnames(w)))
        w = wavname;
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

% Search for m-file containing string wname
wfiltFiles = dir(fullfile(ltfatbasepath,sprintf('wavelets/%s%s.m',wprefix,lower(wavname{1}))));
if(isempty(wfiltFiles))
   error('%s: Unknown wavelet type: %s',upper(mfilename),name); 
else
   % if found, crop '.m' from the filename 
   tmpFile = wfiltFiles.name(1:end-2); 
end

% Synthesis filters delay
d = [];
wfiltNargout = nargout(tmpFile);
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

%END FORMATFILTERS

function w = updateTransDirect(do_ana,w)
if(do_ana)
    w.filts = w.h;
else
    w.filts = w.g;
end


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
        tmph =cellh{ff};
        tmphLen = length(tmph);
        ecent = sum((1:tmphLen-1).*tmph(2:end).^2)/sum(tmph.^2);
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
 
    








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
%   the function produces can (and should) be directly passed to all
%   functions instead of the cell-arrays with wavelet filters which all
%   `wfilt_` functions produces.
%
%   The structure has the following fields:
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
%   See also: fwt, ifwt, waveletfb, multid, wfilt_db

%TO DO: do some chaching
persistent cachw;
persistent cachwDesc;


%Wavelet filters functions definition prefix
wprefix = 'wfilt_';
% Output structure definition
% w = struct('h',{},...
%            'g',{},...
%            'a',[],...
%            'type',[],...
%            'ext',[]);
w.h = {};
w.g = {};
w.a = [];
w.filts = {};

if nargin<1
  return;
end;

% process other parameters

%definput.import = {'fwtcommon','fwt'};
definput.import = {'fwtcommon'};
[flags,kv]=ltfatarghelper({},definput,varargin);

% WILL PROBABLY GO AWAY....
% if(flags.do_ext_null)
%    w.ext = 'per'; 
% else
%    w.ext =  flags.ext; 
% end
%%%%%%%%%%%%%%%%%%%%%%%%%%

if ischar(wavname)
    wname = {wavname};
elseif iscell(wavname)
    wname = wavname;
    if ~ischar(wname{1})
       if iscell(wname{1})
          if(length(wname)==1)
             if(flags.do_ana)
                w.h = wname{1}; 
             else
                w.g = wname{1}; 
             end 
          elseif(length(wname)==2)
             w.h = wname{1};  
             w.g = wname{2};
          else
             error('%s: Unrecognizer format of the filterbank definition.',upper(mfilename)); 
          end
       elseif isnumeric(wname{1})
          if(flags.do_ana)
             w.h = wname; 
          else
             w.g = wname; 
          end 
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

[w.h, w.g, w.a] = feval(tmpFile,wname{2:end});
w.filts = w.h;

% overwrite a if explicitly defined
w.a = formata(length(w.h),kv.a,w.a);
cachw = w;




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




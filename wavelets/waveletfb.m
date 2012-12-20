function [w] = waveletfb(wavname,varargin)
%WAVELETFB  Basic Wavelet Filters
%   Usage:  w = waveletfb(wavname,...);
%
%   Input parameters:
%         wavname : Wavelet filters generating function name (without the prefix).
%
%   Output parameters:
%         w       : Structure defining the filters.
%
%   `w=waveletfb(wavname,...)` produces structure describing one level perfect
%   reconstruction wavelet-type filterbank analysis and synthesis parts.
%   The *wavname* can be string defining function name (without the prefix)
%   to be called or it can be cell-array with the first element beeing the
%   function name and the remaining elements are passed further to the wfilt_ function.
%
%   Function is a wrapper for calling all the functions starting with wfilt_
%   defined in the LTFAT wavelets directory. The structure the function produces can (and should)
%   be directly passed to all functions instead of the cell-arrays with wavelet filters which all wfilt_ functions produces.    
%
%   The structure have the following fields:
%   w.h - analysis filter bank
%   w.g - synthesis filter bank
%   w.a - implicit subsampling factors
%   w.type - dec, undec
%   w.ext - extension type
%
%   See also: fwt, ifwt, waveletfb, multid, wfilt_db

%Wavelet filters functions definition prefix
wprefix = 'wfilt_';

if nargin<1
  error('%s: Too few input parameters.',upper(mfilename));
end;

if ischar(wavname)
    wname = {wavname};
elseif iscell(wavname)
    wname = wavname;
    if ~ischar(wname{1})
    error(['%s: First element of the first agument must be a string denoting the type of ' ...
         'wavelet.'],upper(mfilename));
    end;
else
    error('%s: First argument must be a string or cell.',upper(mfilename));
end;

% Search for m-file containing string wname
wfiltFiles = dir(fullfile(ltfatbasepath,sprintf('wavelets/%s%s*.m',wprefix,lower(wavname{1}))));
if(isempty(wfiltFiles))
   error('%s: Unknown wavelet type: %s',upper(mfilename),name); 
else
   % if found, crop '.m' from the filename 
   tmpFile = wfiltFiles.name(1:end-2); 
end

[w.h, w.g, w.a] = feval(tmpFile,wname{2:end});


% process other parameters
definput.import = {'fwt'};
[flags,kv]=ltfatarghelper({},definput,{varargin{2:end}});

if(flags.do_type_null)
   w.type = 'dec'; 
 else
   w.type = flags.type; 
end
if(flags.do_ext_null)
   w.ext = 'per'; 
else
   w.ext =  flags.ext; 
end


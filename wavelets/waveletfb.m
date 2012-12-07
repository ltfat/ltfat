function [w] = waveletfb(wavname,varargin)
%WAVELETFB   Wavelet Filterbank
%
%
% `w=waveletfb(name,...)` produces structure describing one level perfect
% reconstruction wavelet-type filterbank analysis and synthesis parts.
%
% The structure have the following fields:
% w.h - analysis filter bank
% w.g - synthesis filter bank
% w.a - implicit subsampling factors
% w.type - dec, undec
% w.ext - extension type

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



%ord = 1;
%if(~isempty(varargin)) ord = varargin{1}; end;

% Search for m-file containing string wname
wfiltFiles =dir(fullfile(ltfatbasepath,sprintf('wavelets/%s*.m',wprefix)));
found = 0;
for ii = 1:length(wfiltFiles)
    tmpFile = lower(wfiltFiles(ii).name);
    subtmpFile = tmpFile(length(wprefix)+1:end-2);
    if(strcmp(subtmpFile,wname{1}))
       found = 1; break;
    end
end

if(~found)
   error('%s: Unknown wavelet type: %s',upper(mfilename),name); 
else
   % if found, crop '.m' from the filename 
   tmpFile = tmpFile(1:end-2); 
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



% switch(wname)
%     case('algmband')
%         % ord = 1,2
% 
%        [w.h, w.g, w.a] = wfilt_algmband(ord);
%     case('db')
%        ord = 1;
%        definput.import = {'fwt'};
%        [flags,kv,J]=ltfatarghelper({},definput,{varargin{2:end}});
%        if(~isempty(varargin)) ord = varargin{1}; end;
%        [w.h, w.g, w.a] = wfilt_db(ord);
%        if(flags.do_type_null)
%          w.type = 'dec'; 
%        else
%          w.type = flags.type; 
%        end
%        if(flags.do_ext_null)
%          w.ext = 'per'; 
%        else
%          w.ext =  flags.ext; 
%        end
%     case('dden')
%        no = 1;
%        definput.import = {'fwt'};
%        [flags,kv]=ltfatarghelper({},definput,{varargin{2:end}});
%        if(~isempty(varargin)) no = varargin{1}; end;
%        [w.h, w.g, w.a] = wfilt_dden(no);
%        
%        if(flags.do_type_null)
%          w.type = 'dec'; 
%        else
%          w.type = flags.type; 
%        end
% 
%        if(flags.do_ext_null)
%          w.ext = 'per'; 
%        else
%          w.ext =  flags.ext; 
%        end
%        
%     case('dgrid')
%        no = 1;
%        definput.import = {'fwt'};
%        [flags,kv]=ltfatarghelper({},definput,{varargin{2:end}});
%        if(~isempty(varargin)) no = varargin{1}; end;
%        [w.h, w.g, w.a] = wfilt_dgrid(no);
%        
%        if(flags.do_type_null)
%          w.type = 'dec'; 
%        else
%          w.type = flags.type; 
%        end
% 
%        if(flags.do_ext_null)
%          w.ext = 'per'; 
%        else
%          w.ext =  flags.ext; 
%        end
%        
%      case('dtree')
%        no = 1;
%        definput.import = {'fwt'};
%        [flags,kv]=ltfatarghelper({},definput,{varargin{:}});
%        if(~isempty(varargin)) no = varargin{1}; end;
%        [w.h, w.g, w.a] = wfilt_dtree;
%        
%        w.type = 'dtdwt'; 
% 
% 
%        if(flags.do_ext_null)
%          w.ext = 'per'; 
%        else
%          w.ext =  flags.ext; 
%        end
%        
%      case('hden')
%        no = 1;
%        definput.import = {'fwt'};
%        [flags,kv]=ltfatarghelper({},definput,{varargin{2:end}});
%        if(~isempty(varargin)) no = varargin{1}; end;
%        [w.h, w.g, w.a] = wfilt_hden(no);
%        
%        w.type = 'hddwt'; 
% 
%        if(flags.do_ext_null)
%          w.ext = 'per'; 
%        else
%          w.ext =  flags.ext; 
%        end
%        
%      case('optfs')
%        no = 1;
%        definput.import = {'fwt'};
%        [flags,kv]=ltfatarghelper({},definput,{varargin{2:end}});
%        if(~isempty(varargin)) no = varargin{1}; end;
%        [w.h, w.g, w.a] = wfilt_optfs(no);
%        
%        if(flags.do_type_null)
%          w.type = 'dec'; 
%        else
%          w.type = flags.type; 
%        end
% 
%        if(flags.do_ext_null)
%          w.ext = 'per'; 
%        else
%          w.ext =  flags.ext; 
%        end
%        
%     case('symds')
%        no = 1;
%        definput.import = {'fwt'};
%        [flags,kv]=ltfatarghelper({},definput,{varargin{2:end}});
%        if(~isempty(varargin)) no = varargin{1}; end;
%        [w.h, w.g, w.a] = wfilt_symds(no);
%        
%        if(flags.do_type_null)
%          w.type = 'dec'; 
%        else
%          w.type = flags.type; 
%        end
% 
%        if(flags.do_ext_null)
%          w.ext = 'per'; 
%        else
%          w.ext =  flags.ext; 
%        end
% 
%         
%     otherwise
%         error('%s: Unknown wavelet type: %s',upper(mfilename),name); 
% end;



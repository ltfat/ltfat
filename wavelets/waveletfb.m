function [w, a] = waveletfb(wname,varargin)
%WAVELETFB   Wavelet Filterbank
%
%
% `w=waveletfb(name,...)` produces structure describing one level perfect
% reconstruction wavelet-type filterbank analysis and synthesis parts.
%

if nargin<1
  error('%s: Too few input parameters.',upper(mfilename));
end;

if ~ischar(wname)
  error(['%s: First agument must be a string denoting the type of ' ...
         'wavelet.'],upper(mfilename));
end;



switch(wname)
    case('db')
        ord = 1;
        definput.keyvals.J=[]; 
        definput.import = {'fwt'};
       [flags,kv,J]=ltfatarghelper({'J'},definput,{varargin{2:end}});
       if(isempty(J)) J=1; end
       if(~isempty(varargin)) ord = varargin{1}; end;
       [w.h, w.g, a] = wfilt_db(ord,J);
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
       
    case('dden')
       no = 1;
       definput.import = {'fwt'};
       [flags,kv]=ltfatarghelper({},definput,{varargin{2:end}});
       if(~isempty(varargin)) no = varargin{1}; end;
       [w.h, w.g, a] = wfilt_dden(no);
       
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
       
    case('dgrid')
       no = 1;
       definput.import = {'fwt'};
       [flags,kv]=ltfatarghelper({},definput,{varargin{2:end}});
       if(~isempty(varargin)) no = varargin{1}; end;
       [w.h, w.g, a] = wfilt_dgrid(no);
       
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
       
     case('dtree')
       no = 1;
       definput.import = {'fwt'};
       [flags,kv]=ltfatarghelper({},definput,{varargin{:}});
       if(~isempty(varargin)) no = varargin{1}; end;
       [w.h, w.g, a] = wfilt_dtree;
       
       w.type = 'dtdwt'; 


       if(flags.do_ext_null)
         w.ext = 'per'; 
       else
         w.ext =  flags.ext; 
       end
       
     case('hden')
       no = 1;
       definput.import = {'fwt'};
       [flags,kv]=ltfatarghelper({},definput,{varargin{2:end}});
       if(~isempty(varargin)) no = varargin{1}; end;
       [w.h, w.g, a] = wfilt_hden(no);
       
       w.type = 'hddwt'; 

       if(flags.do_ext_null)
         w.ext = 'per'; 
       else
         w.ext =  flags.ext; 
       end
       
     case('optfs')
       no = 1;
       definput.import = {'fwt'};
       [flags,kv]=ltfatarghelper({},definput,{varargin{2:end}});
       if(~isempty(varargin)) no = varargin{1}; end;
       [w.h, w.g, a] = wfilt_optfs(no);
       
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
       
    case('symds')
       no = 1;
       definput.import = {'fwt'};
       [flags,kv]=ltfatarghelper({},definput,{varargin{2:end}});
       if(~isempty(varargin)) no = varargin{1}; end;
       [w.h, w.g, a] = wfilt_symds(no);
       
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

        
    otherwise
        error('%s: Unknown wavelet type: %s',upper(mfilename),name); 
end;
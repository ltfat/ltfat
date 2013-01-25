function f=iwfbt(c,wt,Ls,varargin)
%WTREE   Inverse Wavelet Filterbank Tree
%
%
% `f=iwtree(c,wtdual)` 

if nargin<2
   error('%s: Too few input parameters.',upper(mfilename));
end;


%% PARSE INPUT
definput.keyvals.Ls=[];    
definput.import = {'fwt'};

if(iscell(c))
    [flags,kv]=ltfatarghelper({},definput,varargin);
    if(flags.do_type_null)
       flags.type = 'dec'; 
    end

    if(flags.do_ext_null)
       flags.ext = 'per'; 
    end
else
    error('%s: Unrecognized coefficient format.',upper(mfilename));
end

f = comp_iwfbt(c,wt,Ls,flags.type,flags.ext);
function f=iwfbt(c,wt,Ls,varargin)
%IWFBT   Inverse Wavelet Filterbank Tree
%
%
%   `f=iwfbt(c,wtdual)` XXX

if nargin<3
   error('%s: Too few input parameters.',upper(mfilename));
end;


%% PARSE INPUT
definput.keyvals.Ls=[];    
definput.import = {'fwt','wfbtcommon'};

if(~iscell(c))
    error('%s: Unrecognized coefficient format.',upper(mfilename));
end

[flags,kv]=ltfatarghelper({},definput,varargin);

% Initialize the wavelet tree structure
wt = wfbtinit(wt,flags.treetype,'syn');

f = comp_iwfbt(c,wt,Ls,'dec',flags.ext);

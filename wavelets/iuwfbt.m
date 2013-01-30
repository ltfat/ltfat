function f=iuwfbt(c,wt,varargin)
%IWFBT   Inverse Undecimated Wavelet Filterbank Tree
%
%
% `f=iwfbt(c,wtdual)` 

if nargin<2
   error('%s: Too few input parameters.',upper(mfilename));
end;


%% PARSE INPUT
definput.keyvals.Ls=[];    
definput.import = {'fwt','wfbtcommon'};

if(~iscell(c))
    error('%s: Unrecognized coefficient format.',upper(mfilename));
end

definput.keyval.Ls = [];
[flags,kv,Ls]=ltfatarghelper({'Ls'},definput,varargin);

% Initialize the wavelet tree structure
wt = wfbtinit(wt,flags.treetype,'syn');


% Estimate output signal length from the number of coefficients
if isempty(Ls)
   [sigHalfLen,W] = size(c{end});
   if(strcmp(flags.ext,'per'))
      Ls = sigHalfLen;  
   else
      error('Not done yet!');
   end
end

f = comp_iwfbt(c,wt,Ls,'undec',flags.ext);
function f=iuwfbt(c,wt,varargin)
%IUWFBT   Inverse Undecimated Wavelet Filterbank Tree
%   Usage:  f = iuwfbt(c,wt,...) 
%
%   `f=iuwfbt(c,wt)` computes XXX 
%
%
%
%

if nargin<2
   error('%s: Too few input parameters.',upper(mfilename));
end;


%% PARSE INPUT
definput.keyvals.Ls=[];    
definput.import = {'wfbtcommon'};


definput.keyval.Ls = [];
[flags,kv,Ls]=ltfatarghelper({'Ls'},definput,varargin);

% Initialize the wavelet tree structure
wt = wfbtinit(wt,flags.treetype,'syn');

%% ----- step 2 : Prepare input parameters
wtPath = nodesBForder(wt,'rev');
nodesUps = nodeFiltUps(wtPath,wt);
rangeLoc = rangeInLocalOutputs(wtPath,wt);
rangeOut = rangeInOutputs(wtPath,wt);
%% ----- step 3 : Run computation
f = comp_iuwfbt(c,wt.nodes(wtPath),nodesUps,rangeLoc,rangeOut);
function f=iwfbt(c,par,varargin)
%IWFBT   Inverse Wavelet Filterbank Tree
%   Usage:  f=iwfbt(c,info);
%           f=iwfbt(c,wt,Ls,...);
%
%   Input parameters:
%         c     : Coefficients stored in a cell-array.
%         wt    : Wavelet Filterbank tree
%         Ls    : Length of the reconstructed signal.
%
%   Output parameters:
%         f     : Reconstructed data.
%
%   `f=iwfbt(c,wt)` reconstructs signal *f* from coefficients *c* using the
%   wavelet filterbank tree *wt*.
%
%   The following flag groups are supported:
%
%         'per','zpd','sym','symw','asym','asymw','ppd','sp0'
%                Type of the boundary handling.
%
%         'dwt','full'
%                Type of the tree to be used.
%
%         'freq','nat'
%                Frequency or natural order of the coefficient subbands.
%
%   Please see the help on |fwt| for a description of the flags.
%
%   Examples:
%   ---------
%   
%   A simple example showing perfect reconstruction using the "full decomposition" wavelet tree:::
% 
%     f = gspi;
%     J = 7;
%     wtdef = {'db10',J,'full'};
%     c = wfbt(f,wtdef);
%     fhat = iwfbt(c,wtdef,length(f));
%     % The following should give (almost) zero
%     norm(f-fhat)
%
%   See also: wfbt, wfbtinit
%

if nargin<2
   error('%s: Too few input parameters.',upper(mfilename));
end;

if(~iscell(c))
   error('%s: Unrecognized coefficient format.',upper(mfilename));
end

if(isstruct(par)&&isfield(par,'fname'))
   if nargin>2
      error('%s: Too many input parameters.',upper(mfilename));
   end
   wt = wfbtinit(par.wt,par.fOrder,'syn');
   Ls = par.Ls;
   ext = par.ext;
   L = wfbtlength(Ls,wt,ext);
else
   if nargin<3
      error('%s: Too few input parameters.',upper(mfilename));
   end;

   %% PARSE INPUT
   definput.keyvals.Ls=[];    
   definput.import = {'fwt','wfbtcommon'};

   [flags,kv,Ls]=ltfatarghelper({'Ls'},definput,varargin);
   ext = flags.ext;
   % Initialize the wavelet tree structure
   wt = wfbtinit(par,flags.forder,'syn');
   % Determine next legal input data length.
   L = wfbtlength(Ls,wt,ext);
end

%% ----- step 3 : Run computation
wtPath = nodesBForder(wt,'rev');
outLengths = nodeInLen(wtPath,L,strcmpi(ext,'per'),wt);
outLengths(end) = Ls;
rangeLoc = rangeInLocalOutputs(wtPath,wt);
rangeOut = rangeInOutputs(wtPath,wt);
f = comp_iwfbt(c,wt.nodes(wtPath),outLengths,rangeLoc,rangeOut,ext);




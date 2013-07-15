function f=iuwfbt(c,par,varargin)
%IUWFBT   Inverse Undecimated Wavelet Filterbank Tree
%   Usage: f = iuwfbt(c,info) 
%          f = iuwfbt(c,wt) 
%
%   Input parameters:
%         c      : Coefficients stored in $L \times M$ matrix.
%         info,w : Transform parameters struct/Wavelet tree definition.
%
%   Output parameters:
%         f     : Reconstructed data.
%
%   `f = iuwfbt(c,info)` reconstructs signal *f* from the coefficients *c* 
%   using parameters from `info` struct. both returned by the |uwfbt| function.
%
%   `f = iuwfbt(c,wt)` reconstructs signal *f* from the wavelet coefficients 
%   *c* using the undecimated wavelet filterbank tree described by *wt* using
%   the "a-trous" algorithm.
%
%   Examples:
%   ---------
%   
%   A simple example showing perfect reconstruction using the "full decomposition" wavelet tree:::
% 
%     f = greasy;
%     J = 6;
%     c = uwfbt(f,{'db8',J,'full'});
%     fhat = iuwfbt(c,{'db8',J,'full'});
%     % The following should give (almost) zero
%     norm(f-fhat)
%
%   See also:  uwfbt, plotwavelets

if nargin<2
   error('%s: Too few input parameters.',upper(mfilename));
end;

if(isstruct(par)&&isfield(par,'fname'))
   if nargin>2
      error('%s: Too many input parameters.',upper(mfilename));
   end
   wt = wfbtinit({'dual',par.wt},par.fOrder);
else
   %% PARSE INPUT
   definput.import = {'wfbtcommon'};
   flags=ltfatarghelper({},definput,varargin);

   % Initialize the wavelet tree structure
   wt = wfbtinit(par,flags.forder);
end


%% ----- step 2 : Prepare input parameters
wtPath = nodesBForder(wt,'rev');
nodesUps = nodeFiltUps(wtPath,wt);
rangeLoc = rangeInLocalOutputs(wtPath,wt);
rangeOut = rangeInOutputs(wtPath,wt);
%% ----- step 3 : Run computation
f = comp_iuwfbt(c,wt.nodes(wtPath),nodesUps,rangeLoc,rangeOut);
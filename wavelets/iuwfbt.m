function f=iuwfbt(c,par,varargin)
%IUWFBT   Inverse Undecimated Wavelet Filterbank Tree
%   Usage: f = iuwfbt(c,info)
%          f = iuwfbt(c,wt)
%
%   Input parameters:
%         c       : Coefficients stored in $L \times M$ matrix.
%         info,wt : Transform parameters struct/Wavelet tree definition.
%
%   Output parameters:
%         f     : Reconstructed data.
%
%   `f = iuwfbt(c,info)` reconstructs signal *f* from the coefficients *c*
%   using parameters from `info` struct. both returned by the |uwfbt| 
%   function.
%
%   `f = iuwfbt(c,wt)` reconstructs signal *f* from the wavelet coefficients
%   *c* using the undecimated wavelet filterbank tree described by *wt*.
%
%   Please see help for |wfbt| description of possible formats of *wt*.
%
%   Filter scaling:
%   ---------------
%
%   As in |uwfbt|, the function recognizes three flags controlling scaling
%   of the filters:
%
%      'sqrt'
%               Each filter is scaled by `1/sqrt(a)`, there *a* is the hop
%               factor associated with it. If the original filterbank is
%               orthonormal, the overall undecimated transform is a tight
%               frame.
%               This is the default.
%
%      'noscale'
%               Uses filters without scaling.
%
%      'scale'
%               Each filter is scaled by `1/a`.
%
%   If 'noscale' is used, 'scale' must have been used in |uwfbt| (and vice
%   versa) in order to obtain a perfect reconstruction.
%
%   Examples:
%   ---------
%
%   A simple example showing perfect reconstruction using the "full
%   decomposition" wavelet tree:::
%
%     f = greasy;
%     J = 6;
%     c = uwfbt(f,{'db8',J,'full'});
%     fhat = iuwfbt(c,{'db8',J,'full'});
%     % The following should give (almost) zero
%     norm(f-fhat)
%
%   See also:  uwfbt, plotwavelets

complainif_notenoughargs(nargin,2,'IUWFBT');

if(isstruct(par)&&isfield(par,'fname'))
   complainif_toomanyargs(nargin,2,'IUWFBT');

   wt = wfbtinit({'dual',par.wt},par.fOrder);
   scaling = par.scaling;

   % Use the "oposite" scaling
   if strcmp(scaling,'scale')
       scaling = 'noscale';
   elseif strcmp(scaling,'noscale')
       scaling = 'scale';
   end
else
   %% PARSE INPUT
   definput.import = {'wfbtcommon','uwfbtcommon'};
   flags=ltfatarghelper({},definput,varargin);

   % Initialize the wavelet tree structure
   wt = wfbtinit(par,flags.forder);
   scaling = flags.scaling;
end


%% ----- step 2 : Prepare input parameters
[nodesBF, rangeLoc, rangeOut] = treeBFranges(wt,'rev');
nodesUps = nodesFiltUps(nodesBF,wt);

%% ----- step 3 : Run computation
f = comp_iuwfbt(c,wt.nodes(nodesBF),nodesUps,rangeLoc,rangeOut,scaling);

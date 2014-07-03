function [g,a] = wfbt2filterbank( wt, varargin)
%WFBT2FILTERBANK  WFBT equivalent non-iterated filterbank
%   Usage: [g,a] = wfbt2filterbank(wt)
%
%   Input parameters:
%         wt : Wavelet filter tree definition
%
%   Output parameters:
%         g   : Cell array containing filters
%         a   : Vector of sub-/upsampling factors
%
%   `[g,a]=wfbt2filterbank(wt)` calculates the impulse responses *g* and the
%   subsampling factors *a* of non-iterated filterbank, which is equivalent
%   to the wavelet filterbank tree described by *wt*. The returned
%   parameters can be used directly in |filterbank| ant other routines.
%
%   `[g,a]=wfbt2filterbank({w,J,'dwt'})` doest the same for the DWT (|FWT|)
%   filterbank tree.
%
%   The function internally calls |wfbtinit| and passes *wt* and all
%   additional parameters to it.
%
%   Examples:
%   ---------
%
%   The following two examples create a multirate identity filterbank
%   using a tree of depth 3. In the first example, the filterbank is
%   identical to the DWT tree:::
%
%      [g,a] = wfbt2filterbank({'db10',3,'dwt'});
%      filterbankfreqz(g,a,1024,'plot','linabs','posfreq');
%
%   In the second example, the filterbank is identical to the full
%   wavelet tree:::
%
%      [g,a] = wfbt2filterbank({'db10',3,'full'});
%      filterbankfreqz(g,a,1024,'plot','linabs','posfreq');
%
%   See also: wfbtinit


complainif_notenoughargs(nargin,1,'WFBT2FILTERBANK');

% build the tree
wt = wfbtinit({'strict',wt},varargin{:});

% Pick just nodes with outputs
wtPath = 1:numel(wt.nodes);
wtPath(nodesOutputsNo(1:numel(wt.nodes),wt)==0)=[];

rangeLoc = nodesLocOutRange(wtPath,wt);
rangeOut = nodesOutRange(wtPath,wt);
[g,a] = nodesMultid(wtPath,rangeLoc,rangeOut,wt);


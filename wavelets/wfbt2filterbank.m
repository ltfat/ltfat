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
%   to the wavelet filterbank tree described by *wt* used in |wfbt|. The 
%   returned parameters can be used directly in |filterbank| and other routines.
%
%   `[g,a]=wfbt2filterbank({w,J,'dwt'})` does the same for the DWT (|FWT|)
%   filterbank tree.
%
%   Please see help on |wfbt| for description of *wt* and help on |fwt| for
%   description of *w* and *J*. 
%
%   The function additionally support the following flags:
%
%   `'freq'`(default),`'nat'`
%     The filters are ordered to produce subbands in the same order as 
%     |wfbt| with the same flag.
%
%   `'scaling_notset'`(default),`'noscale'`,`'scale'`,`'sqrt'`
%     Support for scaling flags as described in |uwfbt|. By default,
%     the returned filterbank *g* and *a* is equivalent to |wfbt|,
%     passing any of the non-default flags results in a filterbank 
%     equivalent to |uwfbt| i.e. scaled and with `a(:)=1`.
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

% AUTHOR: Zdenek Prusa


complainif_notenoughargs(nargin,1,'WFBT2FILTERBANK');

definput.import = {'uwfbtcommon', 'wfbtcommon'};
definput.importdefaults = {'scaling_notset'};
flags = ltfatarghelper({},definput,varargin);

% build the tree
wt = wfbtinit({'strict',wt},flags.forder);

% Pick just nodes with outputs
% wtPath = 1:numel(wt.nodes);
% wtPath(nodesOutputsNo(1:numel(wt.nodes),wt)==0)=[];
% 
% rangeLoc = nodesLocOutRange(wtPath,wt);
% rangeOut = nodesOutRange(wtPath,wt);

[nodesBF, rangeLoc, rangeOut] = treeBFranges(wt);
slice = ~cellfun(@isempty,rangeOut); % Limit to nodes with unconnected outputs
[g,a] = nodesMultid(nodesBF(slice),rangeLoc(slice),rangeOut(slice),wt);

if ~flags.do_scaling_notset
   g = comp_filterbankscale(g,a,flags.scaling);
   a = ones(numel(g),1);
end


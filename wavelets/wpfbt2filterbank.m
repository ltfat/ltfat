function [g,a] = wpfbt2filterbank( wt, varargin)
%WPFBT2FILTERBANK  WPFBT equivalent non-iterated filterbank
%   Usage: [g,a] = wpfbt2filterbank(wt)
%
%   Input parameters:
%         wt : Wavelet filter tree definition
%
%   Output parameters:
%         g   : Cell array containing filters
%         a   : Vector of sub-/upsampling factors
%
%   `wpfbt2filterbank(wt)` calculates the impulse responses *g* and the
%   subsampling factors *a* of non-iterated filterbank, which is equivalent
%   to the wavelet packet filterbank tree described by *wt*. The returned
%   parameters can be used directly in |filterbank|, |ufilterbank| or
%   |filterbank|.
%
%   Please see help on |wfbt| for description of *wt*. The function
%   additionally support the following flags:
%
%   `'freq'`(default),`'nat'`
%      The filters are ordered to produce subbands in the same order as 
%      |wpfbt| with the same flag.
%
%   `'intsqrt'`(default),`'intnoscale'`, `'intscale'`
%      The filters in the filterbank tree are scaled to reflect the
%      behavior of |wpfbt| and |iwpfbt| with the same flags.
%
%   `'scaling_notset'`(default),`'noscale'`,`'scale'`,`'sqrt'`
%     Support for scaling flags as described in |uwpfbt|. By default,
%     the returned filterbank *g* and *a* is equivalent to |wpfbt|,
%     passing any of the non-default flags results in a filterbank 
%     equivalent to |uwpfbt| i.e. scaled and with `a(:)=1`.
%
%   Examples:
%   ---------
%
%   The following two examples create a multirate identity filterbank
%   using a tree of depth 3. In the first example, the filterbank is
%   identical to the DWT tree:::
%
%     [g,a] = wpfbt2filterbank({'db10',3,'dwt'});
%     filterbankfreqz(g,a,1024,'plot','linabs','posfreq');
%
%
%   In the second example, the filterbank is identical to the full
%   wavelet tree:::
%
%     [g,a] = wpfbt2filterbank({'db10',3,'full'});
%     filterbankfreqz(g,a,1024,'plot','linabs','posfreq');
%
%   See also: wfbtinit

% AUTHOR: Zdenek Prusa


complainif_notenoughargs(nargin,1,'WPFBT2FILTERBANK');

definput.import = {'wfbtcommon','uwfbtcommon'};
definput.importdefaults = {'scaling_notset'};
definput.flags.interscaling={'intsqrt','intnoscale','intscale'};
[flags]=ltfatarghelper({},definput,varargin);

% build the tree
wt = wfbtinit({'strict',wt},flags.forder);

wt = comp_wpfbtscale(wt,flags.interscaling);

nIdx = nodesLevelsBForder(wt);
% Now we need to walk the tree by levels
g = {};
a = [];
for ii=1:numel(nIdx)
    rangeLoc = cellfun(@(eEl) 1:numel(eEl.h),wt.nodes(nIdx{ii}),...
                       'UniformOutput',0);
    rangeOut = cellfun(@(eEl) numel(eEl.h),wt.nodes(nIdx{ii}));
    rangeOut = mat2cell(1:sum(rangeOut),1,rangeOut);
    [gtmp,atmp] = nodesMultid(nIdx{ii},rangeLoc,rangeOut,wt);
    g(end+1:end+numel(gtmp)) = gtmp;
    a(end+1:end+numel(atmp)) = atmp;
end
g = g(:);
a = a(:);

if ~flags.do_scaling_notset
   g = comp_filterbankscale(g,a,flags.scaling);
   a = ones(numel(g),1);
end


function nodesIdxs = nodesLevelsBForder(treeStruct)


%find root
nodeNo = find(treeStruct.parents==0);
toGoTrough = [nodeNo];
nodesIdxs = {nodeNo};
inLevel = [1];
counter = 0;
level = 2;
chIdxSum = 0;
while ~isempty(toGoTrough)
   chtmp = find(treeStruct.children{toGoTrough(1)}~=0);
   chIdxtmp = treeStruct.children{toGoTrough(1)}(chtmp);
   counter = counter + 1;

   if(length(nodesIdxs)<level&&~isempty(chIdxtmp))
       nodesIdxs = {nodesIdxs{:},[]}; 
   end
   
   chIdxSum = chIdxSum + length(chIdxtmp);
   if(~isempty(chIdxtmp))
       nodesIdxs{level} = [nodesIdxs{level},chIdxtmp];
   end
   
   toGoTrough = [toGoTrough(2:end),chIdxtmp];

   if(counter==inLevel(level-1))
       counter = 0;
       inLevel(level) = chIdxSum;
       level = level + 1;
       chIdxSum = 0;
   end
end





















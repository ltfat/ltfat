function [nodesBF, rangeLoc, rangeOut] = treeBFranges(wt,varargin)
%TREEBFRANGES Tree nodes output ranges in BF order
%   Usage: [nodesBF, rangeLoc, rangeOut] = treeBFranges(wt);
%          [nodesBF, rangeLoc, rangeOut] = treeBFranges(wt,'rev');
%
%   Input parameters:
%         wt       : Filterbank tree struct.
%   Output parameters:
%         nodesBF  : All nodes in a breadth-first order
%         rangeLoc : Local ranges of unconnected (terminal) outputs
%         rangeOut : Global ranges of unconnected (terminal) outputs
%
%   `[nodesBF, rangeLoc, rangeOut] = treeBFranges(wt)` is a helper function
%   extracting all nodes of a tree in a BF order (root and low-pass first) 
%   (numeric array of indexes `nodesBF`), and two cell arrays of ranges of
%   outputs. Each element of `rangeLoc` specifies range of unconnected
%   outputs of a node with at the corresponding position in `nodesBF`.
%   Elements `rangeOut` specify contain the resulting global indexes 
%   (in the resulting coefficient cell array) of such unconnected nodes.
%
%   `[nodesBF, rangeLoc, rangeOut] = treeBFranges(wt,'rev')` does the same 
%   but the arrays are reversed.


nodesBF = nodeBForder(0,wt);
do_rev = 0;
if ~isempty(varargin(strcmp('rev',varargin)));
   nodesBF = fliplr(nodesBF); 
   do_rev = 1;
end

rangeLoc = nodesLocOutRange(nodesBF,wt);

rangeOut = treeOutRange(wt);
if do_rev
    %rangeOut = nodesOutRange(nodesBF,wt);
    rangeOut = rangeOut(end:-1:1);
end


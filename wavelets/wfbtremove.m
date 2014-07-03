function wt = wfbtremove(d,kk,wt,varargin)
%WFBTREMOVE Remove node(s) from the filterbank tree
%   Usage:  wt = wbftremove(d,kk,wt);
%           wt = wfbtremove(d,kk,wt,'force');
%
%   Input parameters:
%           d   : Level in the tree (0 - root).
%           kk  : Index of the node at level *d* (starting at 0) or array 
%                 of indexes. 
%           wt  : Wavelet filterbank tree structure (as returned from
%                 |wfbtinit|).
%
%   Output parameters:
%           wt : Modified filterbank structure.
%   
%   `wfbtremove(d,kk,wt)` removes existing node at level *d* and index *kk*
%   from the filterbank tree structure *wt*. The function fails if the 
%   node has any children (it is not a leaf node).
%
%   `wfbtremove(d,k,wt,'force')` does the same, but any childern of the
%   node are removed too.
%

complainif_notenoughargs(nargin,3,'WFBTREMOVE');

definput.flags.force = {'noforce','force'};
flags=ltfatarghelper({},definput,varargin);

if isempty(wt.nodes)
   error('%s: Tree is empty.',mfilename); 
end

for k=kk
   [nodeNo,nodeChildIdx] = depthIndex2NodeNo(d,k,wt);
   if(nodeNo==0)
       % removing root 
       rootNo = find(wt.parents==0);
       % check for any children of the root
       if any(wt.children{rootNo}~=0) && ~flags.do_force
            error(['%s: Deleting root node. To delete the whole tree ',...
                   'use FORCE option.'],mfilename,d,k); 
       else
           wt = nodeSubtreeDelete(rootNo,wt);
           continue;
       end
   end

   % check if node exists
   childrenIdx = find(wt.children{nodeNo}~=0);
   found = find(childrenIdx==nodeChildIdx,1);

   if(isempty(found))
        error('%s: Such node (depth=%d, idx=%d) does not exist.',mfilename,d,k); 
   end

   nodeToDelete = wt.children{nodeNo}(nodeChildIdx);
   % Check if it is a leaf (terminal node)
   if any(wt.children{nodeToDelete}~=0) && ~flags.do_force
       error(['%s: Deleting a non-leaf node. To delete whole subtree use ',...
              'FORCE option.'],mfilename);
   else
       wt = nodeSubtreeDelete(nodeToDelete,wt);
   end
end


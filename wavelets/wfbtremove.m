function wt = wfbtremove(d,k,wt,varargin)
%WFBTREMOVE Remove node from the filterbank tree
%   Usage:  wt = wbftremove(d,k,wt);
%           wt = wfbtremove(d,k,wt,'force');
%
%   Input parameters:
%           d   : Level in the tree (0 - root).
%           k   : Index of the node at level *d* (starting at 0).
%           wt  : Wavelet filterbank tree structure (as returned from
%                 |wfbtinit|_).
%
%   Output parameters:
%           wt : Modified filterbank structure.
%   
%   `wfbtremove(d,k,wt)` removes existing node at level *d* and index *k*
%   from the filterbank tree structure *wt*. The function fails if the 
%   node has any children (is not a leaf node).
%
%   `wfbtremove(d,k,wt,'force')` does the same, but any childern of the
%   node are removed too.
%

if(nargin<3)
   error('%s: Too few input parameters.',upper(mfilename)); 
end

definput.flags.force = {'noforce','force'};
[flags,kv]=ltfatarghelper({},definput,varargin);

if(isempty(wt.nodes))
   error('%s: Tree is empty.',mfilename); 
end

[nodeNo,nodeChildIdx] = depthIndex2NodeNo(d,k,wt);
if(nodeNo==0)
    % removing root 
    rootNo = find(wt.parents==0);
    % check for any children of the root
    if(isempty(find(wt.children{rootNo}~=0,1)))
        wt = wtree_init();
    else
        if(flags.do_force)
            wt = wfbtinit();
        else
            error('%s: Deleting root node. To delete the whole tree use FORCE option.',mfilename,d,k); 
        end
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
if(~isempty(find(wt.children{nodeToDelete}~=0,1)))
    if(flags.do_force)
        wt = deleteSubtree(nodeToDelete,wt);
        return;
    else
        error('%s: Deleting non-leaf node. To delete whole subtree use FORCE option.',mfilename);
    end
else
    wt = deleteNode(nodeToDelete,wt); 
end

%wtree = deleteNode(nodeToDelete,wtree);

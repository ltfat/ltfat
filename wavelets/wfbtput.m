function wt = wfbtput(d,k,w,wt,varargin)
%WFBTPUT  Put node to the filterbank tree
%   Usage:  wt = wfbtput(d,k,w,wt);
%           wt = wfbtput(d,k,w,wt,'force');
%
%   Input parameters:
%           d   : Level in the tree (0 - root).
%           k   : Index of the node at level *d* (starting at 0).
%           w   : Node, basic wavelet filterbank.
%           wt  : Wavelet filterbank tree structure (as returned from
%                 |wfbtinit|).
%
%   Output parameters:
%           wt : Modified filterbank structure.
%
%   `wfbtput(d,k,w,wt)` puts the basic filterbank *filt* to the filter
%   tree structure *wt* at level *d* and index *k*. The output is a
%   modified tree structure. *d* and *k* have to specify unconnected output
%   of the leaf node. Error is issued if *d* and *k* points to already
%   existing node. For possible formats of parameter *w* see help of |fwt|.
%   Parameter *wt* has to be a structure returned by |wfbtinit|.
%
%   `wfbtput(d,k,w,wt,'force')` does the same but replaces node at *d* and *k*
%   if it already exists. If the node to be replaced has any children, 
%   the number of outputs of the replacing node have to be equal to number of
%   outputs of the node beeing replaced.
%
  
if(nargin<4)
   error('%s: Too few input parameters.',upper(mfilename)); 
end
definput.flags.force = {'noforce','force'};
[flags,kv]=ltfatarghelper({},definput,varargin);

node = fwtinit(w);

[nodeNo,nodeChildIdx] = depthIndex2NodeNo(d,k,wt);

if(nodeNo==0)
    % adding root 
    if(~isempty(find(wt.parents==0,1)))
        if(flags.do_force)
           rootId = find(wt.parents==0,1);
           % if root has children, check if the new root has the same
           % number of them
           if(~isempty(find(wt.children{rootId}~=0,1)))
              if(length(w)~=length(wt.nodes{rootId}))
                 error('%s: The replacing root have to have %d filters.',mfilename,length(wt.nodes{rootId})); 
              end
           end
        else
            error('%s: Root already defined. Use FORCE option to replace.',mfilename);  
        end
        wt.nodes{rootId} = node;
        return;
    end
    wt.nodes{end+1} = node;
    wt.parents(end+1) = nodeNo;
    wt.children{end+1} = [];
    return;
end

childrenIdx = find(wt.children{nodeNo}~=0);
found = find(childrenIdx==nodeChildIdx,1);
if(~isempty(found))
   if(flags.do_force)
     %check if childrenIdx has any children
     tmpnode = wt.children{nodeNo}(found);  
     if(~isempty(find(wt.children{tmpnode}~=0, 1)))
         if(length(w)~=length(wt.nodes{tmpnode}))
            error('%s: The replacing node have to have %d filters.',mfilename,length(wt.nodes{childrenIdx})); 
         end
     end
     wt.nodes{tmpnode} = node;
     %wtree.a{tmpnode} = a;
     return;
   else
       error('%s: Such node (depth=%d, idx=%d) already exists. Use FORCE option to replace.',mfilename,d,k); 
   end
end


wt.nodes{end+1} = node;
wt.parents(end+1) = nodeNo;
wt.children{end+1} = [];
wt.children{nodeNo}(nodeChildIdx) = numel(wt.parents);
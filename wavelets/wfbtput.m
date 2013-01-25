function wtree = wfbtput(d,k,filt,wtree,varargin)
%WTREE_PUT Put node to the filterbank tree
%   Usage:  wtree = wtreeput(d,k,filt,wtree);
%           wtree = wtreeput(d,k,filt,wtree,'force','a',a);
%
%   Input parameters:
%         d     : 
%         k     :
%         filts :
%         a     :
%         wtree :
%
%   Output parameters:
%         wtree : Modified input structure.
%
%   `wtreeput(d,k,filt,a,wtree)`
%
  
if(nargin<4)
   error('%s: Too few input parameters.',upper(mfilename)); 
end
definput.flags.force = {'noforce','force'};
definput.import = {'fwtcommon'};
[flags,kv]=ltfatarghelper({},definput,varargin);

node = fwtinit(filt,'a',kv.a,flags.ansy);

[nodeNo,nodeChildIdx] = depthIndex2NodeNo(d,k,wtree);

if(nodeNo==0)
    % adding root 
    if(~isempty(find(wtree.parents==0,1)))
        if(flags.do_force)
           rootId = find(wtree.parents==0,1);
           % if root has children, check if the new root has the same
           % number of them
           if(~isempty(find(wtree.children{rootId}~=0,1)))
              if(length(filt)~=length(wtree.nodes{rootId}))
                 error('%s: The replacing root have to have %d filters.',mfilename,length(wtree.nodes{rootId})); 
              end
           end
        else
            error('%s: Root already defined. Use FORCE option to replace.',mfilename);  
        end
        wtree.nodes{rootId} = node;
        return;
    end
    wtree.nodes{end+1} = node;
    wtree.parents(end+1) = nodeNo;
    wtree.children{end+1} = [];
    return;
end

childrenIdx = find(wtree.children{nodeNo}~=0);
found = find(childrenIdx==nodeChildIdx,1);
if(~isempty(found))
   if(doForce)
        %check if childrenIdx has any children
     tmpnode = wtree.children{nodeNo}(found);  
     if(~isempty(find(wtree.children{tmpnode}~=0, 1)))
         if(length(filt)~=length(wtree.nodes{tmpnode}))
            error('%s: The replacing root have to have %d filters.',mfilename,length(wtree.nodes{childrenIdx})); 
         end
     end
     wtree.nodes{tmpnode} = node;
     %wtree.a{tmpnode} = a;
     return;
   else
       error('%s: Such node (depth=%d, idx=%d) already exists. Use FORCE option to replace.',mfilename,d,k); 
   end
end


wtree.nodes{end+1} = node;
wtree.parents(end+1) = nodeNo;
wtree.children{end+1} = [];
wtree.children{nodeNo}(nodeChildIdx) = length(wtree.parents);
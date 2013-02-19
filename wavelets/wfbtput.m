function wfb = wfbtput(d,k,filt,wfb,varargin)
%WFBTPUT  Put node to the filterbank tree
%   Usage:  wtree = wfbtput(d,k,filt,wfb);
%           wtree = wfbtput(d,k,filt,wfb,'force','a',a);
%
%   Output parameters:
%         wtree : Modified input structure.
%
%   `wfbtput(d,k,filt,wtree)` puts the basic filterbank *filt* to the filter
%   tree structure *wtree* at level *d* and index *k*. The output is a
%   modified tree structure.
%
%   `wfbtput(d,k,filt,wfb,'force','a',a)' does the same but XXX Explain this.
%
  
if(nargin<4)
   error('%s: Too few input parameters.',upper(mfilename)); 
end
definput.flags.force = {'noforce','force'};
definput.import = {'fwtcommon'};
[flags,kv]=ltfatarghelper({},definput,varargin);

node = fwtinit(filt,'a',kv.a,flags.ansy);

[nodeNo,nodeChildIdx] = depthIndex2NodeNo(d,k,wfb);

if(nodeNo==0)
    % adding root 
    if(~isempty(find(wfb.parents==0,1)))
        if(flags.do_force)
           rootId = find(wfb.parents==0,1);
           % if root has children, check if the new root has the same
           % number of them
           if(~isempty(find(wfb.children{rootId}~=0,1)))
              if(length(filt)~=length(wfb.nodes{rootId}))
                 error('%s: The replacing root have to have %d filters.',mfilename,length(wfb.nodes{rootId})); 
              end
           end
        else
            error('%s: Root already defined. Use FORCE option to replace.',mfilename);  
        end
        wfb.nodes{rootId} = node;
        return;
    end
    wfb.nodes{end+1} = node;
    wfb.parents(end+1) = nodeNo;
    wfb.children{end+1} = [];
    return;
end

childrenIdx = find(wfb.children{nodeNo}~=0);
found = find(childrenIdx==nodeChildIdx,1);
if(~isempty(found))
   if(doForce)
        %check if childrenIdx has any children
     tmpnode = wfb.children{nodeNo}(found);  
     if(~isempty(find(wfb.children{tmpnode}~=0, 1)))
         if(length(filt)~=length(wfb.nodes{tmpnode}))
            error('%s: The replacing root have to have %d filters.',mfilename,length(wfb.nodes{childrenIdx})); 
         end
     end
     wfb.nodes{tmpnode} = node;
     %wtree.a{tmpnode} = a;
     return;
   else
       error('%s: Such node (depth=%d, idx=%d) already exists. Use FORCE option to replace.',mfilename,d,k); 
   end
end


wfb.nodes{end+1} = node;
wfb.parents(end+1) = nodeNo;
wfb.children{end+1} = [];
wfb.children{nodeNo}(nodeChildIdx) = length(wfb.parents);
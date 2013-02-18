function outRange = rangeInOutputs(nodeNo,treeStruct)
%RANGEINOUTPUTS Index range of the outputs
%   Usage:  outRange = rangeInOutputs(nodeNo,treeStruct);
%
%   Input parameters:
%         nodeNo     : Node index.
%         treeStruct : Structure containing description of the filter tree.
%
%   Output parameters:
%         outRange   : Coefficients stored in J+1 cell-array.
%
%   `rangeInOutputs(nodeNo,treeStruct)` Returns index range in the global
%   tree outputs associated with node `nodeNo`. Empty matrix is returned if
%   node has all outputs connected do children nodes. For definition of the
%   structure see `wfbinit`.
%
%   See also: wfbtinit
%

nodesCount = length(nodeNo);
if(nodesCount>1)
   outRange = cell(nodesCount,1); 
else
    outRange = [];
end

for jj=1:nodesCount
    outRangeTmp = rangeInNodeOutputs(nodeNo(jj),treeStruct);

    if(isempty(outRangeTmp))
        continue;
    end
    %rootId = find(treeStruct.parents==0);
    rootId = nodeNo(jj);
    higherNodes = [];
    while treeStruct.parents(rootId)
         parId = treeStruct.parents(rootId);
          % save idx of all higher nodes
         ch = treeStruct.children{parId};
         childIdx = find(ch==rootId);
         higherNodes(end+1:end+(childIdx-1))=ch(1:childIdx-1);
         rootId = parId;
    end

    noOutPrev = 0;
    for ii=1:length(higherNodes)
        if(higherNodes(ii)==0) 
           noOutPrev=noOutPrev+1; 
        else
           noOutPrev = noOutPrev + noOfSubtreeOutputs(higherNodes(ii),treeStruct);
        end
    end

    outRangeTmp = outRangeTmp + noOutPrev;
    
   if(iscell(outRange))
      outRange{jj} = outRangeTmp;
   else
      outRange = outRangeTmp;
   end
end
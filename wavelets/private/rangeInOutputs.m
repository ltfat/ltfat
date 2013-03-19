function outRange = rangeInOutputs(nodeNo,wt)
%RANGEINOUTPUTS Index range of the outputs
%   Usage:  outRange = rangeInOutputs(nodeNo,treeStruct);
%
%   Input parameters:
%         nodeNo     : Node index.
%         wt         : Structure containing description of the filter tree.
%
%   Output parameters:
%         outRange   : Coefficients stored in J+1 cell-array.
%
%   `rangeInOutputs(nodeNo,wt)` Returns index range in the global
%   tree outputs associated with node `nodeNo`. Empty matrix is returned if
%   node has all outputs connected do children nodes. For definition of the
%   structure see `wfbinit`.
%
%   See also: wfbtinit
%


nodesBFo = nodesBForder(wt,'rev');
nOuts = noOfNodeOutputs(nodesBFo,wt);
nSubtOut = zeros(numel(nodesBFo),1);
for ii=1:numel(nodesBFo)
   childTmp = wt.children{nodesBFo(ii)};
   childTmp(childTmp==0)=[];
   nSubtOut(nodesBFo(ii)) = nOuts(ii) + sum(nSubtOut(childTmp));
end

nodesCount = length(nodeNo);
outRange = cell(nodesCount,1); 

% tic;
% noOut = cellfun(@(nEl) numel(nEl.filts), wt.nodes(nodeNo)) -...
%         cellfun(@(chEl) numel(chEl(chEl~=0)), wt.children(nodeNo));
% toc;

t = 0;
for jj=1:nodesCount
   %tic;
   outRangeTmp = rangeInNodeOutputs(nodeNo(jj),wt);
   %t=t+toc;


    if(isempty(outRangeTmp))
        continue;
    end
    %rootId = find(treeStruct.parents==0);
    rootId = nodeNo(jj);
    higherNodes = [];

%tic;
    while wt.parents(rootId)
         parId = wt.parents(rootId);
          % save idx of all higher nodes
         ch = wt.children{parId};
         childIdx = find(ch==rootId);
         higherNodes(end+1:end+(childIdx-1))=ch(1:childIdx-1);
         rootId = parId;
    end
%t=t+toc;
%tic;
    noOutPrev = 0;
    for ii=1:length(higherNodes)
        if(higherNodes(ii)==0) 
           noOutPrev=noOutPrev+1; 
        else
          % noOutPrev = noOutPrev + noOfSubtreeOutputs(higherNodes(ii),wt);
          noOutPrev = noOutPrev + nSubtOut(higherNodes(ii));
        end
    end
%t=t+toc;

    outRange{jj} = outRangeTmp + noOutPrev;
end

%disp(sprintf('Spend %d ms.',1000*t));



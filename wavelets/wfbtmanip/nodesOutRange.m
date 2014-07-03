function outRange = nodesOutRange(nodeNo,wt)
%NODESOUTRANGE Index range of the outputs
%   Usage:  outRange = nodesOutRange(nodeNo,treeStruct);
%
%   Input parameters:
%         nodeNo     : Node index.
%         wt         : Structure containing description of the filter tree.
%
%   Output parameters:
%         outRange   : Subband idx range.
%
%   `nodesOutRange(nodeNo,wt)` Returns index range in the global
%   tree subbands associated with node `nodeNo`. Empty matrix is returned if
%   node has all outputs connected do children nodes. For definition of the
%   structure see `wfbinit`.
%
%   See also: wfbtinit
%


nodesBFo = fliplr(nodeBForder(0,wt));
nOuts = nodesOutputsNo(nodesBFo,wt);
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

function outRange = rangeInNodeOutputs(nodeNo,treeStruct)
%RANGEINNODEOUTPUTS Index range of the node outputs
%   Usage:  outRange = rangeInNodeOutputs(nodeNo,treeStruct)
%
%   Input parameters:
%         nodeNo     : Node index.
%         treeStruct : Structure containing description of the filter tree.
%
%   Output parameters:
%         outRange   : Index range. 
%
%   `rangeInNodeOutputs(nodeNo,treeStruct)` For definition of the structure
%   see `wfbinit`.
%
%   See also: wfbtinit
%
chIdx = find(treeStruct.children{nodeNo}~=0);
chan = numel(treeStruct.nodes{nodeNo}.g);
outNodes = zeros(chan,1);
outNodes(chIdx) = treeStruct.children{nodeNo}(chIdx);

outRangeStart = 0;
outRange = [];
for ii=1:chan
   if(outNodes(ii)==0)
       outRange(end+1) = outRangeStart+1;
       outRangeStart = outRangeStart+1;
   else
      outRangeStart=outRangeStart + noOfSubtreeOutputs(outNodes(ii),treeStruct);
   end
end

function noOut = noOfSubtreeOutputs(nodeNo,wt)

noChildOut = noOfChildOutputs(nodeNo,wt);
chan = numel(wt.nodes{nodeNo}.g);
child = length(find(wt.children{nodeNo}~=0));
noOut = chan -child + noChildOut;


function noOut = noOfChildOutputs(nodeNo,wt)

noOut = 0;
childrenIdx = find(wt.children{nodeNo}~=0);
children = wt.children{nodeNo}(childrenIdx);
for nn=1:length(children)
   chNodeNo = children(nn);
   chan = numel(wt.nodes{chNodeNo}.g); 
   child = numel(find(wt.children{chNodeNo}~=0));
   noOut = noOut + chan - child;
   noOut = noOut + noOfChildOutputs(chNodeNo,wt);
end







function [g,a] = nodesMultid(wtPath,rangeLoc,rangeOut,wt)
%NODESMULTID Filter tree multirate identity filterbank
%

%clean cache
nodePredecesorsMultId();


% number of outputs of the tree
treeOutputs = sum(cellfun(@(rEl) numel(rEl),rangeOut));

g = cell(treeOutputs,1);
a = zeros(treeOutputs,1);

for ii = 1:numel(wtPath)
   hmi = nodePredecesorsMultId(ii,wt);
   locRange = rangeLoc{ii};
   outRange = rangeOut{ii};
   for jj = 1:length(locRange)
      tmpUpsFac = nodeFiltUps(ii,wt);
      tmpFilt = wt.nodes{ii}.g{locRange(jj)};
      g{outRange(jj)} = struct();
                % 
      g{outRange(jj)}.h = conj(flipud(conv2(hmi,comp_ups(tmpFilt.h(:),tmpUpsFac,1))));
      g{outRange(jj)}.offset = 1-numel(g{outRange(jj)}.h)+nodePredecesorsOrig(-tmpFilt.offset,ii,wt);
   end
   atmp = nodeSub(ii,wt);
   a(outRange) = atmp{1}(locRange);
end
        
        
% clean the cache
nodePredecesorsMultId();


function hmi = nodePredecesorsMultId(nodeNo,treeStruct)
% Build multirate identity of nodes preceeding nodeNo
% chache of the intermediate multirate identities
persistent multIdPre;
% if no paramerer passed, clear the cache
if(nargin==0),  multIdPre = {}; return; end
% in case nodePredecesorsMultId with nodeNo was called before

% if(~isempty(multIdPre))
%   if(length(multIdPre)>=nodeNo&&~isempty(multIdPre{pre(jj)}))
%     hmi = multIdPre{nodeNo};
%   end
% end

startIdx = 1;
hmi = [1];
pre = nodePredecesors(nodeNo,treeStruct);
pre = [nodeNo,pre];
for jj = 1:length(pre)
  if(~isempty(multIdPre))
     if(length(multIdPre)>=pre(jj)&&~isempty(multIdPre{pre(jj)}))
       hmi = multIdPre{pre(jj)};
       startIdx = length(pre)+1 -jj;
       break;
     end
  end
end

pre = pre(end:-1:1);

for ii=startIdx:length(pre)-1
    id = pre(ii);
    hcurr = treeStruct.nodes{id}.g{treeStruct.children{id}==pre(ii+1)}.h(:);
    hcurr = comp_ups(hcurr,nodeFiltUps(id,treeStruct),1);
    hmi = conv2(hmi,hcurr);
end

function predori = nodePredecesorsOrig(baseOrig,nodeNo,treeStruct)
pre = nodePredecesors(nodeNo,treeStruct);
pre = pre(end:-1:1);
if(isempty(pre))
 predori = baseOrig;
 return;
end

pre(end+1) = nodeNo;
predori = baseOrig;
for ii=2:length(pre)
    id = pre(ii);
    predori = nodeFiltUps(id,treeStruct)*baseOrig + predori;
end


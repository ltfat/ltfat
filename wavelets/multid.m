function [out,a,origs] = multid( filts, J, varargin)
%MULTID  Creates equivalent one-level multirate identity filterbank
% filts = cell array
%   - row or collumn orientation -> one level filters, ordered from the lowest
%   central fequency to the highest
%   - cell matrix -> J collumns denote filters at each level, -> no need
%   for J?
%   - dual-tree structures encoded how?
% 
% Deafault downsampling factor for each filter is 2.
% 
% key-value pair 'pack' - actual tree shape (encoded somehow)
% flag dual-tree...


if(nargin<2)
    error('%s: Not enough input arguments',upper(mfilename));
end

definput.keyvals.a = [];
definput.flags.ansy = {'ana','syn'};
definput.flags.treetype = {'dwt','full','cust','wpack'};
[flags,kv,a]=ltfatarghelper({'a'},definput,varargin);  
% for purposes of cleaning the cache
nodePredecesorsMultId();

[fR, fC] = size(filts);

% input is structure defined by waveletfb
if(isstruct(filts))
    a = filts.a;
    if(flags.do_ana)
        filts = filts.h;
    else
        filts = filts.g;
    end
end

if(iscell(filts))
  if(isnumeric(filts{1}))
     if(min([fR, fC])==1)
        fNo = max([fR, fC]);
        for ii = 2:fNo
            if(length(filts{ii})~=length(filts{ii-1}))
                error('%s: Filters are not equal length.',upper(mfilename));
            end
        end
        
        % determine impulse response sample with the "zero" index
        flen = length(filts{1});
           if(flags.do_ana)
              origin = floor( flen/2 );
           else
              origin = floor( flen/2 ) -1; 
           end

        % determine subsampling factors for each filter
        subFac = zeros(fNo,1);
        if(isempty(a))
          subFac(:) = fNo;
        else
           if(length(a)~=fNo)
               error('%s: Number of subsamplers and the number of filters differs.',upper(mfilename));
           end
           subFac = a(:); 
        end
        
        if(flags.do_dwt)
           filtTree = dwtFiltTree(filts,subFac,J);
        elseif(flags.do_full)
           filtTree = fullFiltTree(filts,subFac,J);
        elseif(flags.do_wpack)
            %filtTree = wpFiltTree(filts,subFac,J);
        end
            
        %filtTree = custFiltTree(filts,subFac,J);

        %filtTree = fullFiltTree(filts,subFac,J);
        %filtTree = wpFiltTree(filts,subFac,J);

        origs = cell(noOfOutputs(filtTree),1);
        out = cell(size(origs));
        a = zeros(size(origs));
        
        for ii = 1:length(filtTree.nodes)
            if(noOfNodeOutputs(ii,filtTree)==0), continue; end;
            hmi = nodePredecesorsMultId(ii,filtTree);
            locRange = rangeInLocalOutputs(ii,filtTree);
            outRange = rangeInOutputs(ii,filtTree);
            for jj = 1:length(locRange)
                tmpUpsFac = nodeFiltUps(ii,filtTree);
                out{outRange(jj)} = convolve(hmi,ups(filtTree.nodes{ii}{locRange(jj)},tmpUpsFac,1));
                origs{outRange(jj)} = nodePredecesorsOrig(origin,ii,filtTree);
            end
            atmp = nodeSub(ii,filtTree);
            a(outRange) = atmp(locRange);
        end
        

       out = pad(out,origs);
        
 
     end
  end
  elseif(isstruct(filts))
      
      
      
  else
       error('%s: Filterbank definition should be struct or cell.',upper(mfilename));
end


if(nargout>1)
    % calculate the Multirate identity subsampling factors
%  if(min([fR, fC])==1) 
%     a = multidsubf(J,fNo,subFac);
%  else
%     error('NOT implemented');
%  end
end

% clean the cache
nodePredecesorsMultId();


function a=multidsubf(J,fNo,subFac)
    a = zeros(J*(fNo-1)+1,1);  
  
   for fId=1:fNo-1
       a(end+1-fId)=subFac(end+1-fId);
   end
   
   hsub = subFac(1);
   for jj=2:J
      for fId=1:fNo-1
         a(end+1-(jj-1)*(fNo-1)-fId) = hsub*subFac(fNo+1-fId);
      end  
        hsub= hsub*subFac(1);
   end
   a(1) = hsub;



function h = pad(h,origs)
supp = cell(numel(h),1);

for ii=1:numel(h)
    supp{ii} = [-origs{ii},length(h{ii})-origs{ii}-1];
end


for ii=1:numel(h)
 [s] = supp{ii};
 lS = s(1); hS = s(2);
 if ((lS>0) && (hS>0)) || ((lS<0) && (hS<0) || (abs(lS)==abs(hS)))
     % do nothing, no zero index included in the support
     continue;
 end
 
 if abs(hS) <  abs(lS)   
     h{ii} = [h{ii}(:);zeros(abs(lS)-abs(hS),1);];
 elseif abs(hS) >  abs(lS)
     h{ii} = [zeros(abs(hS)-abs(lS),1);h{ii}(:)]; 
 end

end



function filtTree = wpFiltTree(filts,subFac,J)
  filtTree = fullFiltTree(filts,subFac,J);
  


function filtTree = custFiltTree(filts,subFac,J)
  filtTree = fullFiltTree(filts,subFac,J);
  filtTree = deleteSubtree(5,filtTree);
  filtTree = deleteSubtree(7,filtTree);
  

% create DWT-type depth-J tree (just lowpass branch is decomposed further)
function filtTree = dwtFiltTree(filts,subFac,J)
  for jj=1:J
           filtTree.nodes{jj} = filts;
           filtTree.a{jj} = subFac;
           %filtTree.origins{jj} = subFac;
  end
   for jj=1:J-1
           filtTree.children{jj} = [jj+1];
   end
           filtTree.children{J}  = 0;
           filtTree.parents = [0,1:J-1];
 
% create full depth-J tree
function filtTree = fullFiltTree(filts,subFac,J)
     filtsNo = length(filts);
     filtTree.nodes = {};
     filtTree.a = {};
     filtTree.parents = [];
     filtTree.children = {};
     for jj=1:J
         for ii=1:filtsNo^(jj-1)
           filtTree.nodes{end+1} = filts;
           filtTree.a{end+1} = subFac;
         end
     end
     
     filtTree.parents = zeros(length(filtTree.nodes), 1);
     filtTree.parents(1) = 0;
     filtTree.children{1} = 2:filtsNo+1;
     
      firstfIdx = 2;
      nextfIdx = firstfIdx + filtsNo;
      for jj=2:J
          for ii=1:filtsNo^(jj-1)
            idx = firstfIdx + ii -1;
            filtTree.parents(idx) = ceil((idx-1)/filtsNo);
            filtTree.children{idx} = [((ii-1)*filtsNo+nextfIdx):((ii-1)*filtsNo+nextfIdx +filtsNo-1)];
          end
         firstfIdx = nextfIdx;
         nextfIdx = nextfIdx + filtsNo^(jj);
      end
      
     for jj=1:filtsNo^(J-1)
             filtTree.children{end+1-jj} = [];
     end
     
     filtTree.children = filtTree.children';


% node filter upsampling factor
function upsNo = nodeFiltUps(nodeNo,treeStruct)
tmpNodeNo =nodeNo;
upsNo = 1;
while(treeStruct.parents(tmpNodeNo))
    parentNo = treeStruct.parents(tmpNodeNo);
    upsNo=upsNo*treeStruct.a{parentNo}(find(treeStruct.children{parentNo}==tmpNodeNo));
    tmpNodeNo=parentNo;
end

% node subsampling factor
function subNo = nodeSub(nodeNo,treeStruct)
subNo = nodeFiltUps(nodeNo,treeStruct).*treeStruct.a{nodeNo};

%
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

% Build multirate identity of nodes preceeding nodeNo
function hmi = nodePredecesorsMultId(nodeNo,treeStruct)
% chache of the intermediate multirate identities
persistent multIdPre;
% if no paramerer passed, clear the cache
if(nargin==0),  multIdPre = {}; return; end
% in case nodePredecesorsMultId with nodeNo was called before
if(~isempty(multIdPre))
  if(length(multIdPre)>=nodeNo&&~isempty(multIdPre{pre(jj)}))
    hmi = multIdPre{nodeNo};
  end
end

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
%pre(end+1) = nodeNo;

for ii=startIdx:length(pre)-1
    id = pre(ii);
    hcurr = treeStruct.nodes{id}{find(treeStruct.children{id}==pre(ii+1))};
    hcurr = ups(hcurr,nodeFiltUps(id,treeStruct),1);
    hmi = convolve(hmi,hcurr);
    multIdPre{pre(ii+1)} = hmi;
end






% number of outputs of the tree
function noOut = noOfOutputs(treeStruct)
noOut = 0;
for jj =1:length(treeStruct.nodes)
    chan = length(treeStruct.nodes{jj});
    children = length(find(treeStruct.children{jj}~=0));
    noOut = noOut + chan-children;
end

%Number of output of all children of node
function noOut = noOfChildOutputs(nodeNo,treeStruct)
noOut = 0;
childrenIdx = find(treeStruct.children{nodeNo}~=0);
children = treeStruct.children{nodeNo}(childrenIdx);
for nn=1:length(children)
   chNodeNo = children(nn);
   chan = length(treeStruct.nodes{chNodeNo});
   child = length(find(treeStruct.children{chNodeNo}~=0));
   noOut = noOut + chan -child;
   noOut = noOut + noOfChildOutputs(chNodeNo,treeStruct);
end

%Number of outputs of a subtree starting with nodeNo
function noOut = noOfSubtreeOutputs(nodeNo,treeStruct)
noChildOut = noOfChildOutputs(nodeNo,treeStruct);
chan = length(treeStruct.nodes{nodeNo});
child = length(find(treeStruct.children{nodeNo}~=0));
noOut = chan -child + noChildOut;

%Number of outputs of a subtree starting with nodeNo
function noOut = noOfNodeOutputs(nodeNo,treeStruct)
chan = length(treeStruct.nodes{nodeNo});
child = length(find(treeStruct.children{nodeNo}~=0));
noOut = chan -child;

%Nodes preceeding the nodeNo
function pred = nodePredecesors(nodeNo,treeStruct)
pred = [];
tmpNodeNo = nodeNo;
while treeStruct.parents(tmpNodeNo)~=0
   tmpNodeNo = treeStruct.parents(tmpNodeNo);
   pred(end+1) = tmpNodeNo;
end


% indexes of node outputs
function outRange = rangeInLocalOutputs(nodeNo,treeStruct)
chIdx = find(treeStruct.children{nodeNo}~=0);
chan = length(treeStruct.nodes{nodeNo});
outNodes = zeros(chan,1);
outNodes(chIdx) = treeStruct.children{nodeNo}(chIdx);
outRange = find(outNodes==0);


% indexes of node outputs
function outRange = rangeInNodeOutputs(nodeNo,treeStruct)
chIdx = find(treeStruct.children{nodeNo}~=0);
chan = length(treeStruct.nodes{nodeNo});
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

%Range of idxs in ouptut
function outRange = rangeInOutputs(nodeNo,treeStruct)
outRange = rangeInNodeOutputs(nodeNo,treeStruct);

if(isempty(outRange))
    return;
end
%rootId = find(treeStruct.parents==0);
rootId = nodeNo;
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
 
outRange = outRange + noOutPrev;
%outRange = toGray(outRange);

% can delete just leaves
function treeStruct = deleteNode(nodeNo,treeStruct)
if(~isempty(find(treeStruct.children{nodeNo}~=0)))
    error('Deleting non-leave node!');
end

parId = treeStruct.parents(nodeNo);
toZero = find(treeStruct.children{parId}==nodeNo);
treeStruct.children{parId}(toZero) = 0;

newIdx = 1:length(treeStruct.nodes);
newIdx = newIdx(find(newIdx~=nodeNo));
treeStruct.nodes = {treeStruct.nodes{newIdx}};
treeStruct.a = {treeStruct.a{newIdx}};
%treeStruct.origins = {treeStruct.origins{newIdx}};
treeStruct.parents = treeStruct.parents(newIdx); 
treeStruct.children = {treeStruct.children{newIdx}}';

% and all children and parents with higher idx are lessened
 for ii =1:length(treeStruct.children)
     biggerIdx = find(treeStruct.children{ii}>nodeNo);
     treeStruct.children{ii}(biggerIdx) = treeStruct.children{ii}(biggerIdx)-1;
 end
 biggerIdx = find(treeStruct.parents>nodeNo);
 treeStruct.parents(biggerIdx) = treeStruct.parents(biggerIdx)-1;

% delete whole subtree
function treeStruct = deleteSubtree(nodeNo,treeStruct)
toDelete = nodeSubtreeBF(nodeNo,treeStruct);

for ii = length(toDelete):-1:1
  treeStruct = deleteNode(toDelete(ii),treeStruct); 
  biggerIdx = find(toDelete>toDelete(ii));
  toDelete(biggerIdx) = toDelete(biggerIdx) - 1;
end
treeStruct = deleteNode(nodeNo,treeStruct); 
    
% get tree nodes in the Breadth-first order starting with nodeNo
function subtreeIdx = nodeSubtreeBF(nodeNo,treeStruct)
subtreeIdx = [];

children = treeStruct.children{nodeNo}(find(treeStruct.children{nodeNo}~=0));
subtreeIdx(end+1:end+length(children)) = children;

for ii=1:length(children)
   tmpSbIdx = nodeSubtreeBF(children(ii),treeStruct);
   subtreeIdx(end+1:end+length(tmpSbIdx)) = tmpSbIdx;
end

function y=toGray(x)


y = zeros(size(x));
bitmasks = zeros(32,1);
for ii=0:31
    bitmasks(ii+1) = 2^ii;
end

for ii=1:length(x)
xx=x(ii);
 
ki = zeros(32,1);
xxbin = bitand(bitmasks,xx);

xxbin(find(xxbin>0)) = 1;


for kk=1:length(ki)
  ki(kk) = mod(sum(xxbin(kk:end)),2);
end
y(ii) = bitmasks'*ki;
end








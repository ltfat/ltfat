function [out,a] = wfbtmultid( filts, varargin)
%WFBTMULTID  WFBT equivalent non-iterated filterbank
%   Usage: [out,a,origs] = wfbtmultid(filts)
%
%   Input parameters:
%         filts : Wavelet filter tree definition or basic filters
%                 definition.
%
%   Output parameters:
%         c      : Coefficients stored in a cell-array.
%
%   `multid(filts)` calculates the impulse responses of non-iterated
%   noble multirate-identity wavelet filterbank.
%
%   The following flag groups are supported:
%
%         'dwt','full' : 
%
% 

if(nargin<1)
    error('%s: Not enough input arguments',upper(mfilename));
end

definput.import = {'fwtcommon','wfbtcommon'};
[flags,kv]=ltfatarghelper({},definput,varargin);  
%clean the cache
nodePredecesorsMultId();

filtTree = wfbtinit(filts,varargin{:});
treeOutputs = noOfOutputs(filtTree);
out = cell(treeOutputs,1);
a = zeros(treeOutputs,1);
        
        for ii = 1:length(filtTree.nodes)
            if(noOfNodeOutputs(ii,filtTree)==0), continue; end;
            hmi = nodePredecesorsMultId(ii,filtTree);
            locRange = rangeInLocalOutputs(ii,filtTree);
            outRange = rangeInOutputs(ii,filtTree);
            for jj = 1:length(locRange)
                tmpUpsFac = nodeFiltUps(ii,filtTree);
                tmpFilt = filtTree.nodes{ii}.filts{locRange(jj)};
                out{outRange(jj)} = wfiltstruct('FIR');
                out{outRange(jj)}.h = convolve(hmi,ups(tmpFilt.h,tmpUpsFac,1));
                out{outRange(jj)}.d = nodePredecesorsOrig(tmpFilt.d,ii,filtTree);
            end
            atmp = nodeSub(ii,filtTree);
            a(outRange) = atmp(locRange);
        end


% clean the cache
nodePredecesorsMultId();

%no padding used
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
%pre(end+1) = nodeNo;

for ii=startIdx:length(pre)-1
    id = pre(ii);
    hcurr = treeStruct.nodes{id}.filts{find(treeStruct.children{id}==pre(ii+1))}.h;
%     if(doSyn)
%        hcurr = treeStruct.nodes{id}.g{find(treeStruct.children{id}==pre(ii+1))};
%     else
%        hcurr = treeStruct.nodes{id}.h{find(treeStruct.children{id}==pre(ii+1))};
%     end
    hcurr = ups(hcurr,nodeFiltUps(id,treeStruct),1);
    hmi = convolve(hmi,hcurr);
   % multIdPre{pre(ii+1)} = hmi;
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













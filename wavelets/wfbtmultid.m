function [wfb,a] = wfbtmultid( wtdef, varargin)
%WFBTMULTID  WFBT equivalent non-iterated filterbank
%   Usage: [wfb,a] = wfbtmultid(wtdef)
%
%   Input parameters:
%         wtdef : Wavelet filter tree definition or basic filters
%                 definition.
%
%   Output parameters:
%         wfb   : Cell array containing impulse responses
%         a     : Array of sub-/upsampling factors
%
%   `wfbtmultid(wtdef)` calculates the impulse responses of non-iterated
%   noble multirate-identity wavelet filterbank.
%
%   The function internally calls |wfbtinit| and passes all parameters to it.
%   
%   Examples:
%   ---------
%   
%   The following two examples create a multirate identity filterbank
%   using a tree of depth 3. In the first example, the filterbank is
%   identical to the DWT tree:::
%
%     [wfb,a] = wfbtmultid({'db10',3,'dwt'});
%     wtfftfreqz(wfb);
%
%   In the second example, the filterbank is identical to the full
%   wavelet tree:::
%
%    [wfb,a] = wfbtmultid({'db10',3,'full'});
%    wtfftfreqz(wfb);
%
%   See also: wfbtinit


if(nargin<1)
    error('%s: Not enough input arguments',upper(mfilename));
end

definput.import = {'fwtcommon','wfbtcommon'};
[flags,kv]=ltfatarghelper({},definput,varargin);  

%clean cache
nodePredecesorsMultId();

% build the tree
filtTree = wfbtinit(wtdef,varargin{:});

% number of outputs of the tree
treeOutputs = noOfOutputs(filtTree);
wfb = cell(treeOutputs,1);
a = zeros(treeOutputs,1);
        
        for ii = 1:length(filtTree.nodes)
            % skip nodes with no outputs
            if(noOfNodeOutputs(ii,filtTree)==0), continue; end;
            hmi = nodePredecesorsMultId(ii,filtTree);
            locRange = rangeInLocalOutputs(ii,filtTree);
            outRange = rangeInOutputs(ii,filtTree);
            for jj = 1:length(locRange{1})
                tmpUpsFac = nodeFiltUps(ii,filtTree);
                tmpFilt = filtTree.nodes{ii}.g{locRange{1}(jj)};
                wfb{outRange{1}(jj)} = struct('fc',0,'realonly',0);
                % 
                wfb{outRange{1}(jj)}.h = flipud(conv2(hmi,comp_ups(tmpFilt.h(:),tmpUpsFac,1)));
                wfb{outRange{1}(jj)}.offset = -nodePredecesorsOrig(tmpFilt.d,ii,filtTree);
            end
            atmp = nodeSub(ii,filtTree);
            a(outRange{1}) = atmp{1}(locRange{1});
        end

% clean the cache
nodePredecesorsMultId();
% END WFBTMULTID

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













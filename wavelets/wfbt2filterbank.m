function [g,a] = wfbt2filterbank( wtdef, varargin)
%WFBTMULTID  WFBT equivalent non-iterated filterbank
%   Usage: [g,a] = wfbt2filterbank(wtdef)
%
%   Input parameters:
%         wtdef : Wavelet filter tree definition
%
%   Output parameters:
%         g   : Cell array containing filters
%         a   : Vector of sub-/upsampling factors
%
%   `wfbtmultid(wtdef)` calculates the impulse responses *g* and the 
%   subsampling factors *a* of non-iterated filterbank, which is equivalent
%   to the wavelet filterbank tree described by *wtdef*. The returned 
%   parameters can be used directly in |filterbank|, |ufilterbank| or 
%   |filterbank| . 
%   
%   The filters are scaled if *a* is not returned. 
%
%   The function internally calls |wfbtinit| and passes *wtdef* and all 
%   additional parameters to it.   
%   
%   Examples:
%   --------- 
%   
%   The following two examples create a multirate identity filterbank
%   using a tree of depth 3. In the first example, the filterbank is
%   identical to the DWT tree:::
%
%     [g,a] = wfbt2filterbank({'db10',3,'dwt'});
%     
%
%   In the second example, the filterbank is identical to the full
%   wavelet tree:::
%
%     [g,a] = wfbt2filterbank({'db10',3,'full'});
%
%   See also: wfbtinit


if(nargin<1)
    error('%s: Not enough input arguments',upper(mfilename));
end

%clean cache
nodePredecesorsMultId();

% build the tree
filtTree = wfbtinit({'strict',wtdef},varargin{:});

% number of outputs of the tree
treeOutputs = noOfOutputs(filtTree);
g = cell(treeOutputs,1);
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
                g{outRange{1}(jj)} = struct();
                % 
                g{outRange{1}(jj)}.h = conj(flipud(conv2(hmi,comp_ups(tmpFilt.h(:),tmpUpsFac,1))));
                g{outRange{1}(jj)}.offset = 1-numel(g{outRange{1}(jj)}.h)+nodePredecesorsOrig(-tmpFilt.offset,ii,filtTree);
            end
            atmp = nodeSub(ii,filtTree);
            a(outRange{1}) = atmp{1}(locRange{1});
        end

if nargout<2
   % Scale filters if a is not returned
   for nn=1:numel(g)
       g{nn}.h = g{nn}.h/sqrt(a(nn));
   end
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













function [out,a,origs] = wfbtmultid( filts, varargin)
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

    % determine impulse response sample with the "zero" index
    % TODO: do it better
%    if(flags.do_ana)
       flen = length(filtTree.nodes{1}.h{1});
       origin = floor( flen/2 );
%     else
%        flen = length(filtTree.nodes{1}.g{1});
%        origin = floor( flen/2 ) -1; 
%     end


        origs = cell(noOfOutputs(filtTree),1);
        out = cell(size(origs));
        a = zeros(size(origs));
        
        for ii = 1:length(filtTree.nodes)
            if(noOfNodeOutputs(ii,filtTree)==0), continue; end;
            hmi = nodePredecesorsMultId(ii,~flags.do_ana,filtTree);
            locRange = rangeInLocalOutputs(ii,filtTree);
            outRange = rangeInOutputs(ii,filtTree);
            for jj = 1:length(locRange)
                tmpUpsFac = nodeFiltUps(ii,filtTree);
               % tmpFilt = filtTree.nodes{ii}.filts{locRange(jj)};
                tmpFilt = filtTree.nodes{ii}.h{locRange(jj)};
 %               if flags.do_ana
                   out{outRange(jj)} = convolve(hmi,ups(tmpFilt,tmpUpsFac,1));
 %               else
 %                  out{outRange(jj)} = convolve(hmi,ups(filtTree.nodes{ii}.g{locRange(jj)},tmpUpsFac,1)); 
 %               end
                origs{outRange(jj)} = nodePredecesorsOrig(origin,ii,filtTree);
                %origs{outRange(jj)} = nodePredecesorsOrig(tmpFilt.d,ii,filtTree);
            end
            atmp = nodeSub(ii,filtTree);
            a(outRange) = atmp(locRange);
        end
       out = pad(out,origs);

      
      



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


% function a=multidsubf(J,fNo,subFac)
%     a = zeros(J*(fNo-1)+1,1);  
%   
%    for fId=1:fNo-1
%        a(end+1-fId)=subFac(end+1-fId);
%    end
%    
%    hsub = subFac(1);
%    for jj=2:J
%       for fId=1:fNo-1
%          a(end+1-(jj-1)*(fNo-1)-fId) = hsub*subFac(fNo+1-fId);
%       end  
%         hsub= hsub*subFac(1);
%    end
%    a(1) = hsub;



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

function hmi = nodePredecesorsMultId(nodeNo,doSyn,treeStruct)
% Build multirate identity of nodes preceeding nodeNo
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
    if(doSyn)
       hcurr = treeStruct.nodes{id}.g{find(treeStruct.children{id}==pre(ii+1))};
    else
       hcurr = treeStruct.nodes{id}.h{find(treeStruct.children{id}==pre(ii+1))};
    end
    hcurr = ups(hcurr,nodeFiltUps(id,treeStruct),1);
    hmi = convolve(hmi,hcurr);
    multIdPre{pre(ii+1)} = hmi;
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













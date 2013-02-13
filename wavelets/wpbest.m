function [c,wt] = wpbest(f,w,J,varargin)
%WPBEST
%
%

if nargin<3
  error('%s: Too few input parameters.',upper(mfilename));
end;

if ~isnumeric(J) || ~isscalar(J)
  error('%s: "J" must be a scalar.',upper(mfilename));
end;

if(J<1 || rem(J,1)~=0)
   error('%s: J must be a positive integer.',upper(mfilename)); 
end

w = fwtinit(w,'ana');

definput.import = {'fwt'};
definput.flags.buildOrder = {'bottomup','topdown'};
definput.flags.bestWhat = {'tree','level'};
definput.keyvals.entropy = {'shannon'};
[flags,kv]=ltfatarghelper({},definput,varargin);



if(~iscell(kv.entropy))
   kv.entropy = {kv.entropy}; 
end

% test if the chosen entropy measure is additive
do_additive = isAdditive(kv.entropy);

if(flags.do_bottomup)
   % do full-tree Wavelet Packet decomposition beforehand and prune
   wt = wfbtinit({w,J},'full','nat');
   c = wpfbt(f,wt,flags.ext);
   % calculate entropy of all subbands
   cEnt = zeros(length(c),1);
   for ii=1:length(c)
       cEnt(ii) = wentwrap(c{ii},kv.entropy{1},kv.entropy{2:end});
   end
    
   % Nodes in the reverse BF order
   treePath = nodesBForder(wt,'rev');
   [pOutIdxs,chOutIdxs] = rangeWpBF(wt,'rev');
   % Nodes to be removed
   removeNodes = [];
   % omit root
   for ii=1:length(pOutIdxs)-1
       pEnt = cEnt(pOutIdxs(ii));
       chEnt = cEnt(chOutIdxs{ii});
       
       % if parent entropy is smaller than combidend children entropy
       if(pEnt<=sum(chEnt))
           removeNodes(end+1) = treePath(ii);
       else
           if(do_additive)
              % Works just for additive measures
              cEnt(pOutIdxs(ii)) = sum(chEnt);
           else
              % Set parent entropy to value obtanied by concatenating child
              % subbands
              cEnt(pOutIdxs(ii)) = wentwrap({c{chOutIdxs{ii}}},kv.entropy{1},kv.entropy{2:end}); 
           end
       end
   end
   for ii=1:length(removeNodes)
       wt = deleteSubtree(removeNodes(ii),wt);
   end
     
    
elseif(flags.do_topdown)
    error('%s: Flag %s not supported yet.',upper(mfilename),flags.buildOrder);
    % Build the from the root
    wt = wfbtinit({w,1});
    

end

% finally do the analysis using the created best tree
c = wfbt(f,wt,flags.ext,'freq'); 
%END WPBEST

function ad = isAdditive(entFunc)
x = 1:5;
x = x./norm(x);
ent1 = wentwrap(x,entFunc{1},entFunc{2:end});

ent2 = 0;
for ii=1:length(x)
    ent2 = ent2 + wentwrap(x(ii),entFunc{1},entFunc{2:end});
end

ad = ent1==ent2;




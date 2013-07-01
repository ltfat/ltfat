function [c,info] = wpbest(f,w,J,varargin)
%WPBEST  Best Tree selection
%   Usage: c = wpbest(f,w,J,...);
%
%   Input parameters:
%         f   : Input data.
%         w   : Wavelet Filterbank.
%         J   : Maximum depth of the tree.
%
%   Output parameters:
%         c   : Coefficients stored in a cell-array.
%         wt  : Structure defining the best tree.
%
%   `[c,wt]=wpbest(f,w,J)` selects the best sub-tree *wt* from the full
%   tree with max. depth *J*, which minimises the coefficient entropy.
%   First, the depth *J* wavelet packet decomposition is performed. Then
%   the nodes are traversed in the reverse breadth-first and entropy of
%   the node input and combined entropy of the node outputs is compared.
%   If the node input entropy is less than the combined output entropy,
%   the current node and all possible descendant nodes are marked to be 
%   deleted, if not, the input is assigned the combidend output entropy.
%   At the end, the marked nodes are removed and the resulting tree is
%   considered to be a best bases in the chosen entropy sense.
%   
%  
%   Examples:
%   ---------
%   
%   A simple example of calling |wpbest| :::
%     
%     f = gspi;
%     J = 8;
%     [c,info] = wpbest(f,'sym10',J,'entropy',{'wlpnorm',1.3});
%
%     % Use 2/3 of the space for the first plot, 1/3 for the second.
%     subplot(3,3,[1 2 4 5 7 8]);
%     plotwavelets(c,info,44100,90);
%
%     subplot(3,3,[3 6 9]);
%     N=cellfun(@numel,c); L=sum(N); a=L./N;
%     plot(a,'o');
%     xlim([1,numel(N)]);
%     view(90,-90);
%     xlabel('Channel no.');
%     ylabel('Subsampling rate / samples');
%
%   References: wick92lecton tas94near-bestbasis   

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
   % Do full-tree Wavelet Packet decomposition beforehand and prune.
   wt = wfbtinit({w,J,'full'},'nat');
   c = wpfbt(f,wt,flags.ext);
   % calculate entropy of all subbands
   cEnt = zeros(length(c),1);
   for ii=1:length(c)
       cEnt(ii) = wentwrap(c{ii},kv.entropy{1},kv.entropy{2:end});
   end
    
   % Nodes in the reverse BF order
   treePath = nodesBForder(wt,'rev');
   % Relationships between nodes
   [pOutIdxs,chOutIdxs] = rangeWpBF(wt,'rev');
   % Nodes to be removed
   removeNodes = [];
   % Skip root.
   for ii=1:length(pOutIdxs)-1
       pEnt = cEnt(pOutIdxs(ii));
       chEnt = cEnt(chOutIdxs{ii});
       
       % if parent entropy is smaller than combidend children entropy
       if(pEnt<=sum(chEnt))
           removeNodes(end+1) = treePath(ii);
       else
           if(do_additive)
              % Set parent entropy to the sum of the children entropy. 
              cEnt(pOutIdxs(ii)) = sum(chEnt);
           else
              % Set parent entropy to value obtanied by concatenating child
              % subbands.
              cEnt(pOutIdxs(ii)) = wentwrap({c{chOutIdxs{ii}}},kv.entropy{1},kv.entropy{2:end}); 
           end
       end
   end
   % Do tree prunning.
   for ii=1:length(removeNodes)
       wt = deleteSubtree(removeNodes(ii),wt);
   end
     
    
elseif(flags.do_topdown)
    error('%s: Flag %s not supported yet.',upper(mfilename),flags.buildOrder);
    % Build the from the root
    wt = wfbtinit({w,1});
    

end

% finally do the analysis using the created best tree
[c,info] = wfbt(f,wt,flags.ext,'freq'); 
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




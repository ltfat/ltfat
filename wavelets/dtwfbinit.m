function [dtw,info] = dtwfbinit(dtwdef,varargin)
%DTWFBINIT Dual-Tree Wavelet Filterbank Initialization
%   Usage:  dtw=dtwfbinit(dtwdef);
%
%   Input parameters:
%         dtwdef : Dual-tree filterbank definition.
%
%   Output parameters:
%         dtw    : Dual-tree filtarbank structure.
%
%   `dtwfinit()` (a call without aguments) creates an empty structure. It
%   has the same fields as the struct. returned from |wfbtinit| plus a field
%   to hold nodes from the second tree:
%
%      .nodes      Filterbank nodes of the first tree
% 
%      .dualnodes  Filterbank nodes of the second tree
%
%      .children   Indexes of children nodes
%
%      .parents    Indexes of a parent node
%
%      .forder     Frequency ordering of the resultant frequency bands.   
%
%   `dtwfinit({dtw,J,flag})` creates a structure representing a dual-tree
%   wavelet filterbank of depth *J*, using dual-tree wavelet filters 
%   specified by `dtw`. The shape of the tree is controlled by `flag`.
%
%   
%   Additional optional flags
%   -------------------------
%
%   `'freq'`,`'nat'`
%      Frequency or natural ordering of the resulting coefficient subbands.
%      When working with wavelet packet filterbank trees, the subbands are
%      
%      This does not affect a 
%      Default ordering is `'freq'`.



% Output structure definition.
% The structure has the same fields as returned by wfbtinit
% but contains additional field .dualnodes containing
% filters of the dual tree
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dtw = wfbtinit();
dtw.dualnodes = {};
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
info.istight = 0;

if nargin < 1
    return;
end

% Frequency or natural ordering
definput.import = {'wfbtcommon'};
[flags,kv]=ltfatarghelper({},definput,varargin);

% Strict throws an error if filterbank is not tight and ana: or syn:
% is not specified.
do_strict = 0;
% Dual returns the 'other' filterbank. 
do_dual = 0;

% Check 'strict'
if iscell(dtwdef) && ischar(dtwdef{1}) && strcmpi(dtwdef{1},'strict')
   do_strict = 1;
   dtwdef = dtwdef{2:end};
end
% Check 'dual'
if iscell(dtwdef) && ischar(dtwdef{1}) && strcmpi(dtwdef{1},'dual')
   do_dual = 1;
   dtwdef = dtwdef{2:end};
end

% FIRST, check if dtwdef is already a struct
if isstruct(dtwdef)
   if isfield(dtwdef,'nodes') && isfield(dtwdef,'dualnodes')
      dtw = dtwdef;
   
       if do_dual || do_strict
          nodesArg = dtw.nodes;
          
          % Undo the frequency ordering
          if dtw.freqOrder
             dtw = nat2freqOrder(dtw,'rev');
          end
          
          doDualTreeFilt = cellfun(@(nEl) strcmp(nEl.wprefix,'wfiltdt_'),...
                                   nodesArg);
        
          if do_dual
             nodesArg = cellfun(@(nEl) {'dual',nEl},nodesArg,...
                                'UniformOutput',0);
          end
          if do_strict
             nodesArg = cellfun(@(nEl) {'strict',nEl},nodesArg,...
                                'UniformOutput',0);
             
          end
          info.istight = 1;
          dtw.nodes = {};
          dtw.dualnodes = {};
      
         
          rangeRest = 1:numel(nodesArg);
          rangeRest(dtw.parents==0) = []; 
                    
          for ii=rangeRest
              if doDualTreeFilt(ii)
                % This is dual-tree specific filterbank
                [dtw.nodes{ii},infotmp] = fwtinit(nodesArg{ii},'wfiltdt_');
                dtw.dualnodes{ii} = dtw.nodes{ii};
                dtw.nodes{ii}.h(:,2) = [];
                dtw.nodes{ii}.g(:,2) = [];
                dtw.dualnodes{ii}.h(:,1) = [];
                dtw.dualnodes{ii}.g(:,1) = [];
              else
                [dtw.nodes{ii},infotmp] = fwtinit(nodesArg{ii}); 
                dtw.dualnodes{ii} = dtw.nodes{ii};
              end
              info.istight = info.istight && infotmp.istight;
          end
          
          % Treat root separately
          [rootNode,infotmp] = fwtinit(nodesArg{dtw.parents==0}); 
          dtw = replaceRoots(dtw,rootNode);
          info.istight = info.istight && infotmp.istight;
          
          % Do the filter frequency shuffling again, since the filters were
          % overwritten in fwtinit.
          if dtw.freqOrder
             dtw = nat2freqOrder(dtw);
          end
       end

       % Do filter shuffling if flags.do_freq differs from the wt.freqOrder.
       % Frequency and natural oreding coincide for DWT.
       if dtw.freqOrder && ~flags.do_freq
          dtw = nat2freqOrder(dtw,'rev');
          dtw.freqOrder = ~dtw.freqOrder;
       elseif ~dtw.freqOrder && flags.do_freq
          dtw = nat2freqOrder(dtw);
          dtw.freqOrder = ~dtw.freqOrder;
       end   
   else
       error('%s: Invalid dual-tree structure format.',upper(mfilename));
   end
   return;
end

% Parse the other params
% Tree type
definput.flags.treetype = {'dwt','full','doubleband','quadband',...
                           'octaband','root'};
% First stage filterbank
definput.keyvals.first = [];
% Leaf filterbank
definput.keyvals.leaf = [];
% Depth of the tree
definput.keyvals.J = [];
wdef = dtwdef{1};
[flags2,kv2,J]=ltfatarghelper({'J'},definput,dtwdef(2:end));
complainif_notposint(J,'J');

% Now dtwdef is this {dtw,J,flag,'first',w}
if do_dual
   wdef = {'dual',wdef};
end
if do_strict
   wdef = {'strict',wdef};
end

if ~(ischar(wdef) || iscell(wdef))
    error('%s: Unrecognized format of dual-tree filters.',upper(mfilename));
end

% Get the dual-tree filters
[w, dtinfo] = fwtinit(wdef,'wfiltdt_');

info.istight = dtinfo.istight;
info.dw = w;

% Determine the first-stage wavelet filters
if ~isfield(dtinfo,'defaultfirst') && isempty(kv2.first)
    error('%s: No first stage wavelet filters specified.',upper(mfilename));
end


if ~isempty(kv2.first) 
   if do_dual
      kv2.first = {'dual',kv2.first};
   end
   if do_strict
      kv2.first = {'strict',kv2.first};
   end

   [kv2.first, firstinfo] = fwtinit(kv2.first);
   isfirsttight = firstinfo.istight;
else
   kv2.first = dtinfo.defaultfirst;
   isfirsttight = dtinfo.defaultfirstinfo.istight;
end



isleaftight = [];
if ~(flags2.do_dwt || flags2.do_root)
    % Determine leaf nodes (only valid for wavelet packets)
    if ~isfield(dtinfo,'defaultleaf') && isempty(kv2.leaf) 
        error('%s: No leaf wavelet filters specified.',...
              upper(mfilename));
    else
       if isempty(kv2.leaf)
          kv2.leaf = dtinfo.defaultleaf; 
          isleaftight = dtinfo.defaultleafinfo.istight;
       else
          if do_dual
             kv2.leaf = {'dual',kv2.leaf};
          end
          if do_strict
             kv2.leaf = {'strict',kv2.leaf};
          end
          [kv2.leaf, leafinfo] = fwtinit(kv2.leaf);
          isleaftight = leafinfo.istight;
       end
    end
end


% Extract filters for dual trees
% This is a bit clumsy...
w1 = w;
w1.h = w1.h(:,1);
w1.g = w1.g(:,1);

w2 = w;
w2.h = w2.h(:,2);
w2.g = w2.g(:,2);


% Initialize both trees
dtw  = wfbtinit({w1,J,flags2.treetype}, 'nat');
dtw2 = wfbtinit({w2,J,flags2.treetype}, 'nat');
% Merge tree definitions to a single struct.
dtw.dualnodes = dtw2.nodes;
dtw = replaceRoots(dtw,kv2.first);


% Replace the 'packet leaf nodes' (see Bayram)
if ~(flags2.do_dwt || flags2.do_root)
    filtNo = numel(w1.g);
    if flags2.do_doubleband 
       for jj=1:J-1
         dtw = wfbtput(2*(jj+1)-1,1:filtNo-1,kv2.leaf,dtw,'force');
       end
    elseif flags2.do_quadband
       idx = 1:filtNo-1;
       idx = [idx,idx+filtNo];
       dtw = wfbtput(2,idx,kv2.leaf,dtw,'force');
       for jj=1:J-1
         dtw = wfbtput(3*(jj+1)-2,1:filtNo-1,kv2.leaf,dtw,'force');
         dtw = wfbtput(3*(jj+1)-1,1:2*filtNo-1,kv2.leaf,dtw,'force');
       end 
    elseif flags2.do_octaband
       idx = 1:filtNo-1;idx = [idx,idx+filtNo];
       dtw = wfbtput(2,idx,kv2.leaf,dtw,'force');
       idx = 1:2*filtNo-1;idx = [idx,idx+2*filtNo];
       dtw = wfbtput(3,idx,kv2.leaf,dtw,'force');
       for jj=1:J-1
         dtw = wfbtput(4*(jj+1)-3,1:filtNo-1,kv2.leaf,dtw,'force');
         dtw = wfbtput(4*(jj+1)-2,1:2*filtNo-1,kv2.leaf,dtw,'force');
         dtw = wfbtput(4*(jj+1)-1,1:4*filtNo-1,kv2.leaf,dtw,'force');
       end 
    elseif flags2.do_full
       for jj=2:J-1
          idx = 1:filtNo^jj-1;
          idx(filtNo^(jj-1))=[];
          dtw = wfbtput(jj,idx,kv2.leaf,dtw,'force');
       end 
    else
        error('%s: Something is seriously wrong!',upper(mfilename));
    end 
end


% Do filter shuffling if frequency ordering is required,
dtw.freqOrder = flags.do_freq;
if flags.do_freq
   dtw = nat2freqOrder(dtw); 
end

% 
info.istight = isfirsttight && info.istight;
if ~isempty(isleaftight)
   info.istight = info.istight && isleaftight;
end


function dtw = replaceRoots(dtw,rootNode)
% Replace the root nodes

firstTmp = rootNode;
firstTmp.h = cellfun(@(hEl) setfield(hEl,'offset',hEl.offset+1),...
                        firstTmp.h,'UniformOutput',0);
firstTmp.g = cellfun(@(gEl) setfield(gEl,'offset',gEl.offset+1),...
                        firstTmp.g,'UniformOutput',0);
                 
% First tree root 
dtw.nodes{dtw.parents==0} = rootNode;  
% Second tree root (shifted by 1 sample)
dtw.dualnodes{dtw.parents==0} = firstTmp;  




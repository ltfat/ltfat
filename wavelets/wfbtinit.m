function wt = wfbtinit(wtdef,varargin)
%WFBTINIT Initialize Filterbank Tree
%   Usage:  wt = wfbtinit(wtdef);
%
%   Input parameters:
%         wtdef : Filterbank tree definition.
%
%   Output parameters:
%         wt    : Structure describing the filter tree.
%
%   `wfbtinit()` creates empty structure. The structure describing the 
%   tree has the following fields:
%
%     .nodes     Actual impulse responses
%   
%     .children  Indexes of children nodes
% 
%     .parents   Indexes of a parent node
%
%     .forder    Frequency ordering of the resultant frequency bands.
%
%   `wfbt=wfbtinit({w,J,flag})` creates a filterbank tree of depth *J*. The
%   parameter *w* defines a basic wavelet filterbank. For all possible formats
%   see |fwt|.  The following optional flags (still inside of the
%   cell-array) are recognized:
%
%   'dwt','full'
%     Type of the tree to be created.
%
%   Additional flags:
%
%   'freq','nat'
%     Frequency or natural ordering of the coefficient subbands. The direct
%     usage of the wavelet tree (`'nat'` option) does not produce coefficient
%     subbans ordered according to the frequency. To achieve that, some 
%     filter shuffling has to be done (`'freq'` option).  
%
%   See also: wfbtput, wfbtremove

% TO DO: Do some caching

% Output structure definition.
% Effectively, it describes a ADT tree.
% .nodes, .children, .parents ale all arrays of the same length and the j-th
% node in the tree is desribed by wt.nodes{j}, wt.children{j} and
% wt.parents(j). wt.nodes{j} is the actual data stored in the tree node 
% (a structure returned form fwtinit) and wt.parents(j) and wt.children{j}
% define relationship to other nodes. wt.parents(j) is an index of the
% parent in the arrays. wt.children{j} is an array of indexes of the
% children nodes.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
wt.nodes = {}; 
wt.children = {};
wt.parents = [];
wt.freqOrder = 0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% return empty struct if no argument was passed
if(nargin<1)
    return;
end


do_strict = 0;
do_dual = 0;

% Check 'strict'
if iscell(wtdef) && ischar(wtdef{1}) && strcmpi(wtdef{1},'strict')
   do_strict = 1;
   wtdef = wtdef{2:end};
end
% Check 'dual'
if iscell(wtdef) && ischar(wtdef{1}) && strcmpi(wtdef{1},'dual')
   do_dual = 1;
   wtdef = wtdef{2:end};
end

definput.import = {'wfbtcommon'};
[flags,kv]=ltfatarghelper({},definput,varargin);

% If wtdef is already a structure
if (isstruct(wtdef)&&isfield(wtdef,'nodes'))
    wt = wtdef;
    
    if do_dual || do_strict
       nodesArg = wt.nodes;
       if do_dual
          nodesArg = cellfun(@(nEl) {'dual',nEl},nodesArg,'UniformOutput',0);
       end
       if do_strict
          nodesArg = cellfun(@(nEl) {'strict',nEl},nodesArg,'UniformOutput',0);
       end
       wt.nodes = cellfun(@(nEl) fwtinit(nEl),nodesArg,'UniformOutput',0);
       % Do the filter frequency shuffling again, since the filters were
       % overwritten in fwtinit.
       if wt.freqOrder
          wt = nat2freqOrder(wt); 
       end
    end

    % Do filter shuffling if flags.do_freq differs from the wt.freqOrder.
    % Frequency and natural oreding coincide for DWT.
    if wt.freqOrder ~= flags.do_freq
       wt = nat2freqOrder(wt); 
       wt.freqOrder = ~wt.freqOrder;
    end
    return;
end

% break if the input parameter is not in the correct format
if ~(iscell(wtdef)) || isempty(wtdef)
    error('%s: Unsupported filterbank tree definition.',upper(mfilename));
end

% Creating new tree
% Now wtdef is this {w,J,flag}
wdef = wtdef{1};
definput = [];
definput.flags.treetype = {'full','dwt','root'};
definput.keyvals.J = [];
[flags2,kv2,J]=ltfatarghelper({'J'},definput,wtdef(2:end));

if do_dual
   wdef = {'dual',wdef};
end

if do_strict
   wdef = {'strict',wdef};
end

w = fwtinit(wdef);

% Doing one-node tree
if flags2.do_root
   J = 1;
end

if isempty(J)
   error('%s: Unspecified J.',upper(mfilename));
end

if flags2.do_dwt || J==1
   % fill the structure to represent a DWT tree
   for jj=0:J-1
      wt = wfbtput(jj,0,w,wt);
   end
elseif flags2.do_full
   % fill the structure to represent a full wavelet tree
   for jj=0:J-1
      for ii=0:numel(w.g)^(jj)-1
         wt = wfbtput(jj,ii,w,wt);
      end
   end
end

% Do filter shuffling if frequency ordering is required,
if flags.do_freq
   wt = nat2freqOrder(wt); 
   wt.freqOrder = 1;
else
   wt.freqOrder = 0;
end





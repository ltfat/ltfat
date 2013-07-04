function wt = wfbtinit(wtdef,varargin)
%WFBTINIT Initialize Filterbank Tree
%   Usage:  wt = wfbtinit();
%           wt = wfbtinit(wtdef,...);
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
%   `wtree.nodes`
%      actual impulse responses
%   
%   `wtree.children`
%      indexes of children nodes
% 
%   `wtree.parents`
%      indexes of a parent node
%
%   `wfbt=wfbtinit({w,J,flag})` creates filterbank tree of depth *J*. Parameter *w* 
%   defines basic wavelet filterbank. For all possible formats see |fwt|.
%   The following optional flags (still inside of the cell-array) are
%   recognized:
%
%   'dwt','full'
%     Type of the tree to be created.
%
%   The following additional flag groups are supported:
%
%   'freq','nat'
%     Frequency or natural order of the coefficient subbands.
%

% TO DO: Do some caching

% output structure definition
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
wt.nodes = {};
wt.children = {};
wt.parents = [];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% return empty struct if no argument was passed
if(nargin==0)
    return;
end


definput.import = {'wfbtcommon','fwtcommon'};
% Contents of the arg_fwtcommon:
% definput.flags.ansy = {'ana','syn'};
% definput.keyvals.a = [];

% Contents of the arg_wfbtcommon:
% definput.flags.treetype = {'full','dwt'};
% definput.flags.forder = {'freq','nat'};
[flags,kv]=ltfatarghelper({},definput,varargin);

% if filts is a structure describing filterbank tree
if((isstruct(wtdef)&&isfield(wtdef,'nodes')))
    wt = wtdef;
    % set analysis or synthesis filters as the active ones
    for ii=1:numel(wt.nodes)
        if(flags.do_ana)
           wt.nodes{ii}.filts = wt.nodes{ii}.h;
        else
           wt.nodes{ii}.filts = wt.nodes{ii}.g;
        end
    end
    % Do filter shuffling if frequency ordering is required,
    % otherwise, the outputs of the tree will be in the natural order.
    % Frequency and natural oreding coincide for DWT.
    if(flags.do_freq)
       wt = nat2freqOrder(wt); 
    end
    % return modified input structure
    return;
end

% break if the input parameter is not in the correct format
if ~(iscell(wtdef)&&numel(wtdef)>=2&&isnumeric(wtdef{2}))
    error('%s: Unsupported filterbank tree definition.',upper(mfilename));
end


J = wtdef{2};

if(numel(wtdef)>=3)
   varargin2 = wtdef(3);
else
   varargin2 = [];
end

definput = [];
definput.flags.treetype = {'full','dwt'};
flags2=ltfatarghelper({},definput,varargin2);

filtsStruct = fwtinit(wtdef{1},flags.ansy);
filtsNo = numel(filtsStruct.filts);

if(flags2.do_dwt || J==1)
   % fill the structure to represent DWT tree
   for jj=1:J
      wt.nodes{jj} = filtsStruct;
   end
   for jj=1:J-1
      wt.children{jj} = [jj+1];
   end
      wt.children{J}  = 0;
      wt.parents = [0,1:J-1]; 
elseif(flags2.do_full)
   % fill the structure to represent full wavelet tree
   for jj=1:J
      for ii=1:filtsNo^(jj-1)
         wt.nodes{end+1} = filtsStruct;
      end
   end
     
   wt.parents = zeros(filtsNo, 1);
   wt.parents(1) = 0;
   wt.children{1} = 2:filtsNo+1;
     
   firstfIdx = 2;
   nextfIdx = firstfIdx + filtsNo;
   for jj=2:J
      for ii=1:filtsNo^(jj-1)
         idx = firstfIdx + ii -1;
         wt.parents(idx) = ceil((idx-1)/filtsNo);
         wt.children{idx} = [((ii-1)*filtsNo+nextfIdx):((ii-1)*filtsNo+nextfIdx +filtsNo-1)];
      end
      firstfIdx = nextfIdx;
      nextfIdx = nextfIdx + filtsNo^(jj);
   end
      
   for jj=1:filtsNo^(J-1)
      wt.children{end+1-jj} = [];
   end
     
   % I messed something here
   wt.children = wt.children(:)'; 
   wt.parents = wt.parents(:)';
     
else
    %just root node
     wt.nodes{1} = filtsStruct;
     wt.children{1} = zeros(filtsNo,1);
     wt.parents(1) = 0;
end

% Do filter shuffling if frequency ordering is required,
if(flags.do_freq)
   wt = nat2freqOrder(wt); 
end





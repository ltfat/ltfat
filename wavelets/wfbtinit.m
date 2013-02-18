function wfbt = wfbtinit(filts,varargin);
%WFBTINIT Initialize Filterbank Tree
%   Usage:  wfbt = wfbtinit();
%           wfbt = wfbtinit(filts);
%
%   Input parameters:
%         filts   : Basic wavelet filterbank.
%
%   Output parameters:
%         wtree   : Structure describing the filter tree.
%
%   `wfbtinit()` creates empty structure. The structure describing the 
%   tree has the following fields:
%
%     `wtree.nodes`
%        actual impulse responses
%     
%     `wtree.children`
%        indexes of children nodes
% 
%     `wtree.parents`
%        indexes of a parent node
%
%   `wfbt=wfbtinit({w,J})` creates filterbank tree of depth *J*. Parameter `w` 
%   can be either structure obtained form the |fwtinit|_ function or a cell
%   array, whose first element is a name of the function defining the basic
%   wavelet filters (`wfilt_` prefix) and the other elements are parameters
%   passed to the function e.g. : `wfbtinit({{'db',10},4})`.
%
%   The following flag groups are supported:
%
%         'dwt','full'
%                Type of the tree to be created.
%
%         'freq','nat'
%                Frequency or natural order of the coefficient subbands.
%

% TO DO: Do some caching

% output structure definition
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
wfbt.nodes = {};
wfbt.children = {};
wfbt.parents = [];
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
if((isstruct(filts)&&isfield(filts,'nodes')))
    wfbt = filts;
    % set analysis or synthesis filters as the active ones
    nodes = wfbt.nodes;
    for ii=1:length(nodes)
        if(flags.do_ana)
           wfbt.nodes{ii}.filts = wfbt.nodes{ii}.h;
        else
           wfbt.nodes{ii}.filts = wfbt.nodes{ii}.g;
        end
    end
    % Do filter shuffling if frequency ordering is required,
    % otherwise, the outputs of the tree will be in the natural order.
    % Frequency and natural oreding coincide for DWT.
    if(flags.do_freq)
       wfbt = nat2freqOrder(wfbt); 
    end
    % return modified input structure
    return;
end

% if no J was specified use J=1; e.g. function was called like so: wfbtinit({'db',10})
if(iscell(filts)&&ischar(filts{1}))
    filts = {filts,1};
end

% break if the input parameter is not in the correct format
if ~(iscell(filts)&&length(filts)==2&&isnumeric(filts{2}))
    error('%s: Unsupported filterbank tree definition.',upper(mfilename));
end


J = filts{2};
filtsStruct = fwtinit(filts{1},flags.ansy);
filtsNo = length(filtsStruct.filts);

if(flags.do_dwt || J==1)
   % fill the structure to represent DWT tree
   for jj=1:J
      wfbt.nodes{jj} = filtsStruct;
   end
   for jj=1:J-1
      wfbt.children{jj} = [jj+1];
   end
      wfbt.children{J}  = 0;
      wfbt.parents = [0,1:J-1]; 
elseif(flags.do_full)
   % fill the structure to represent full wavelet tree
   for jj=1:J
      for ii=1:filtsNo^(jj-1)
         wfbt.nodes{end+1} = filtsStruct;
      end
   end
     
   wfbt.parents = zeros(filtsNo, 1);
   wfbt.parents(1) = 0;
   wfbt.children{1} = 2:filtsNo+1;
     
   firstfIdx = 2;
   nextfIdx = firstfIdx + filtsNo;
   for jj=2:J
      for ii=1:filtsNo^(jj-1)
         idx = firstfIdx + ii -1;
         wfbt.parents(idx) = ceil((idx-1)/filtsNo);
         wfbt.children{idx} = [((ii-1)*filtsNo+nextfIdx):((ii-1)*filtsNo+nextfIdx +filtsNo-1)];
      end
      firstfIdx = nextfIdx;
      nextfIdx = nextfIdx + filtsNo^(jj);
   end
      
   for jj=1:filtsNo^(J-1)
      wfbt.children{end+1-jj} = [];
   end
     
   wfbt.children = wfbt.children'; 
     
else
    %just root node
     wfbt.nodes{1} = filtsStruct;
     wfbt.children{1} = zeros(filtsNo,1);
     wfbt.parents(1) = 0;
end

% Do filter shuffling if frequency ordering is required,
if(flags.do_freq)
   wfbt = nat2freqOrder(wfbt); 
end





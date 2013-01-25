function wtree = wfbtinit(filts,varargin);
%WFBTINIT Initialize Filterbank Tree
%   Usage:  wtree = wtreeinit();
%           wtree = wtreeinit(filts);
%
%   Input parameters:
%         filts   : Basic wavelet filterbank.
%
%   Output parameters:
%         wtree   : Structure describing the filter tree.
%
%   `wtreeinit()` creates empty structure. The structure describing the 
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
%   `wtreeinit({w,J})` creates filterbank tree of depth *J*. Parameter `w` 
%   can be either structure obtained form the |fwtinit|_ function or a cell
%   array, whose first element is a name of the function defining the basic
%   wavelet filters (|wfilt_|_ prefix) and the other elements are parameters
%   passed to the function e.g. : `wtreeinit({{'db',10},4})`.
%

% TO DO: Do some chaching

% initial empty structure
wtree.nodes = {};
%wtree.a = {};
wtree.children = {};
wtree.parents = [];
%wtree.origins = {};

if(nargin==0)
    return;
end

if((isstruct(filts)&&isfield(filts,'nodes')))
    wtree = filts;
    return;
end

if ~(iscell(filts)&&length(filts)==2&&isnumeric(filts{2}))
    error('%s: Unsupported filterbank tree definition.',upper(mfilename));
end

definput.import = {'wfbtcommon'};

[flags,kv]=ltfatarghelper({},definput,varargin);
% explicitly defined subsampling factors are preferred

J = filts{2};
filtsStruct = fwtinit(filts{1});

filtsNo = max([length(filtsStruct.h),length(filtsStruct.g)]);

if(flags.do_dwt || J==1)
  for jj=1:J
           wtree.nodes{jj} = filtsStruct;
           %wtree.a{jj} = a;
           %filtTree.origins{jj} = subFac;
   end
   for jj=1:J-1
           wtree.children{jj} = [jj+1];
   end
           wtree.children{J}  = 0;
           wtree.parents = [0,1:J-1]; 
           
elseif(flags.do_full)
     for jj=1:J
         for ii=1:filtsNo^(jj-1)
           wtree.nodes{end+1} = filtsStruct;
           %wtree.a{end+1} = a;
         end
     end
     
     wtree.parents = zeros(filtsNo, 1);
     wtree.parents(1) = 0;
     wtree.children{1} = 2:filtsNo+1;
     
      firstfIdx = 2;
      nextfIdx = firstfIdx + filtsNo;
      for jj=2:J
          for ii=1:filtsNo^(jj-1)
            idx = firstfIdx + ii -1;
            wtree.parents(idx) = ceil((idx-1)/filtsNo);
            wtree.children{idx} = [((ii-1)*filtsNo+nextfIdx):((ii-1)*filtsNo+nextfIdx +filtsNo-1)];
          end
         firstfIdx = nextfIdx;
         nextfIdx = nextfIdx + filtsNo^(jj);
      end
      
     for jj=1:filtsNo^(J-1)
             wtree.children{end+1-jj} = [];
     end
     
     wtree.children = wtree.children'; 
     
     % TO DO: Do it better!
     wtree = nat2freqOrder(0,wtree);
     wtree = nat2freqOrder(1,wtree);
else
    %just root node
     wtree.nodes{1} = filtsStruct;
     %wtree.a{1} = a;
     wtree.children{1} = zeros(filtsNo,1);
     wtree.parents(1) = 0;
end





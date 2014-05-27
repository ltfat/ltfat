function [dtw,info] = dtwfbinit(dtwdef,varargin)



% Output structure definition.
% The structure has the same fields as returned by wfbtinit
% but contains additional field .dualnodes containing
% filters of the dual tree
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dtw.nodes = {};
dtw.dualnodes = {};
dtw.children = {};
dtw.parents = [];
dtw.freqOrder = 0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
info.istight = 0;


definput.import = {'wfbtcommon'};
[flags,kv]=ltfatarghelper({},definput,varargin);

% FIRST, check if dtwdef is already a struct or a cell array with one of
% the elements being a struct
if iscell(dtwdef) && ... 
   any(cellfun(@(dEl) isstruct(dEl) && isfield(dtwdef,'dualnodes'),dtwdef)) || ...
   isstruct(dtwdef) && ...
   isfield(dtwdef,'nodes') && isfield(dtwdef,'dualnodes')

   [dtw,info1] = wfbtinit(dtwdef,flags.do_freq);
   if isstruct(dtwdef)
      dtwdef.nodes = dtwdef.dualnodes;
   else
      dtwdwfstruct = dtwdef{isstruct(dtwdef)};
      dtwdef{isstruct(dtwdef)}.nodes = dtwdwfstruct.dualnodes; 
   end
   [dtw2,info2] = wfbtinit(dtwdef,flags.do_freq); 

   dtw.dualnodes = dtw2.nodes;
   info.istight = info1.istight && info2.istight;
   return;
end


do_strict = 0;
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

wdef = dtwdef{1};

definput.flags.treetype = {'dwt','full','doubleband','quadband',...
                           'octaband','root'};
definput.keyvals.first = [];
definput.keyvals.leaf = [];
definput.keyvals.J = [];
definput.flags.swapping = {'noswap','swap'};
[flags2,kv2,J]=ltfatarghelper({'J'},definput,dtwdef(2:end));
complainif_notposint(J,'J');

if flags2.do_swap && ~(flags2.do_dwt)
    error('%s: Filter swapping can be used only with dwt tree.',...
          upper(mfilename));
end

% Now dtwdef is this {dtw,J,flag,'first',w}
% Isolate 

if do_dual
   wdef = {'dual',wdef};
end
if do_strict
   wdef = {'strict',wdef};
end

[w, info] = fwtinit(wdef,'wfiltdt_');

firstinfo = info;
if ~isempty(kv2.first)
   if do_dual
      kv2.first = {'dual',kv2.first};
   end
   if do_strict
      kv2.first = {'strict',kv2.first};
   end
   [kv2.first, firstinfo] = fwtinit(kv2.first);
else
%   if abs(sum(1/w.a) - 1) < 1e-6
      % Use symorth3 
      [kv2.first, firstinfo] = fwtinit('db10');
%   end
end

% This is only relevant for packets
leafinfo = info;
if ~isempty(kv2.leaf) && ~(flags2.do_dwt || flags2.do_root)
   if do_dual
      kv2.leaf = {'dual',kv2.leaf};
   end
   if do_strict
      kv2.leaf = {'strict',kv2.leaf};
   end
   [kv2.leaf, leafinfo] = fwtinit(kv2.leaf);
else
   [kv2.leaf, leafinfo] = fwtinit('db10');
end

% Extract dual trees
w1 = w;
w1.h(end/2+1:end) = [];
w1.g(end/2+1:end) = [];
w1.a(end/2+1:end) = [];

w2 = w;
w2.h(1:end/2) = [];
w2.g(1:end/2) = [];
w2.a(1:end/2) = [];

dtw = wfbtinit({w1,J,flags2.treetype},'nat');

% Replace the root node
 if ~isempty(kv2.first)
     firstTmp = kv2.first;
     %firstTmp.h=cellfun(@(hEl) setfield(hEl,'h',[0;0;hEl.h]),firstTmp.h,'UniformOutput',0);
     firstTmp.h{1}.h=flipud([firstTmp.h{1}.h]);
     firstTmp.h{2}.h=flipud([firstTmp.h{2}.h]);
     
     %firstTmp.g=cellfun(@(gEl) setfield(gEl,'h',[0;gEl.h]),firstTmp.g,'UniformOutput',0);
     firstTmp.g{1}.h=flipud([0;0;firstTmp.g{1}.h]);
     firstTmp.g{2}.h=flipud([0;0;firstTmp.g{2}.h]);
     
     dtw = wfbtput(0,0,firstTmp,dtw,'force');
 end

dtw2 = wfbtinit({w2,J,flags2.treetype},'nat');

 if ~isempty(kv2.first)
     firstTmp = kv2.first;
     %firstTmp.h=cellfun(@(hEl) setfield(hEl,'h',[0;hEl.h]),firstTmp.h,'UniformOutput',0);
     firstTmp.h{1}.h=flipud([0;firstTmp.h{1}.h]);
     firstTmp.h{2}.h=flipud([0;firstTmp.h{2}.h]);
     %firstTmp.g=cellfun(@(gEl) setfield(gEl,'h',[0;gEl.h]),firstTmp.g,'UniformOutput',0);
     firstTmp.g{1}.h=flipud([0;firstTmp.g{1}.h]);
     firstTmp.g{2}.h=flipud([0;firstTmp.g{2}.h]);
     dtw2 = wfbtput(0,0,firstTmp,dtw2,'force');
 end

dtw.dualnodes = dtw2.nodes;
%dtw.dualnodes(1) = dtw.nodes(1);

% Replace the 'packet leaf nodes' (see Bayram)
if ~(flags2.do_dwt || flags2.do_root) && J>2
    if flags2.do_doubleband || flags2.do_quadband || flags2.do_octaband
        % Do something different
    else
       filtNo = numel(w1.h);  
       for jj=2:J-1
          idx = 1:2^jj-1;
          idx(2^(jj-1))=[];
          dtw = wfbtput(jj,idx,kv2.leaf,dtw,'force');
       end 
    end 
end


if flags2.do_swap
    % Here we are working with DWT tree only
    % We assume the nodes in dtw.nodes and dtw.dualnodes
    % to be ordered from root to highest level
    if J>2
       [dtw.nodes(3:end), dtw.dualnodes(3:end)] = ...
           deal(dtw.dualnodes(3:end), dtw.nodes(3:end));
    end
end


% Do filter shuffling if frequency ordering is required,
dtw.freqOrder = flags.do_freq;
if flags.do_freq
   dtw = nat2freqOrder(dtw); 
end

info.istight = firstinfo.istight && info.istight;









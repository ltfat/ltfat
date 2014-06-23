function [dtw,info] = dtwfbinit(dtwdef,varargin)



% Output structure definition.
% The structure has the same fields as returned by wfbtinit
% but contains additional field .dualnodes containing
% filters of the dual tree
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dtw = wfbtinit();
dtw.dualnodes = {};
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

% Determine the first-stage wavelet filters
if ~isfield(info,'defaultfirst') && isempty(kv2.first)
    error('%s: No first stage wavelet filters specified.',upper(mfilename));
end

if isempty(kv2.first)
    kv2.first = info.defaultfirst;
end

if do_dual
   kv2.first = {'dual',kv2.first};
end
if do_strict
   kv2.first = {'strict',kv2.first};
end

[kv2.first, firstinfo] = fwtinit(kv2.first);

leafinfo = [];
if ~(flags2.do_dwt || flags2.do_root)
    % Determine leaf nodes (only valid for wavelet packets)
    if ~isfield(info,'defaultleaf') && isempty(kv2.leaf) 
        error('%s: No leaf wavelet filters specified.',...
              upper(mfilename));
    else
       if isempty(kv2.leaf)
          kv2.leaf = info.defaultleaf; 
       end
       if do_dual
          kv2.leaf = {'dual',kv2.leaf};
       end
       if do_strict
          kv2.leaf = {'strict',kv2.leaf};
       end
       [kv2.leaf, leafinfo] = fwtinit(kv2.leaf);
    end
end


% Extract dual trees
% This is a bit clumsy...
w1 = w;
w1.h(end/2+1:end) = [];
w1.g(end/2+1:end) = [];
w1.a(end/2+1:end) = [];

w2 = w;
w2.h(1:end/2) = [];
w2.g(1:end/2) = [];
w2.a(1:end/2) = [];


% Initialize both trees
dtw  = wfbtinit({w1,J,flags2.treetype}, 'nat');
dtw2 = wfbtinit({w2,J,flags2.treetype}, 'nat');


% Replace the root nodes
if ~isempty(kv2.first)
   firstTmp = kv2.first;
   firstTmp.h = cellfun(@(hEl) setfield(hEl,'offset',hEl.offset+1),...
                        firstTmp.h,'UniformOutput',0);
   firstTmp.g = cellfun(@(gEl) setfield(gEl,'offset',gEl.offset+1),...
                        firstTmp.g,'UniformOutput',0);
   dtw = wfbtput(0,0,firstTmp,dtw,'force');
   dtw2 = wfbtput(0,0,kv2.first,dtw2,'force');
end

dtw.dualnodes = dtw2.nodes;

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
       [dtw.nodes(3:2:end), dtw.dualnodes(3:2:end)] = ...
           deal(dtw.dualnodes(3:2:end), dtw.nodes(3:2:end));
    end
end


% Do filter shuffling if frequency ordering is required,
dtw.freqOrder = flags.do_freq;
if flags.do_freq
   dtw = nat2freqOrder(dtw); 
end

info.istight = firstinfo.istight && info.istight && ...
               ~isempty(leafinfo) && leafinfo.istight;









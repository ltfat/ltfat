function [g,a,info] = dtwfb2filterbank( dualwt, varargin)
%DTWFB2FILTERBANK DTWFB equivalent non-iterated filterbank
%   Usage: [g,a] = dtwfb2filterbank(dualwt)
%          [g,a,info] = dtwfb2filterbank(...)
%
%   Input parameters:
%      dualwt : Dual-tree wavelet filterbank specification.
%
%   Output parameters:
%      g      : Cell array of filters.
%      a      : Downsampling rate for each channel.
%      info   : Additional information.
%
%   `[g,a] = dtwfb2filterbank(dualwt)` constructs a set of filters *g* and
%   subsampling factors *a* of a non-iterated filterbank, which is 
%   equivalent to the dual-tree wavelet filterbank defined by *dualwt*. 
%   The returned parameters can be used directly in |filterbank| and other
%   routines. The format of *dualwt* is the same as in |dtwfb| and 
%   |dtwfbreal|.  
%
%   The function internally calls |dtwfbinit| and passes *dualwt* and all
%   additional parameters to it.
%
%   `[g,a,info] = dtwfb2filterbank(...)` additionally outputs a *info*
%   struct containing equivalent filterbanks of individual real-valued
%   trees as fields *info.g1* and *info.g2*.
%
%   Additional parameters:
%   ----------------------
%
%   `'real'`
%      By default, the function returns a filtebank equivalent to |dtwfb|.
%      The filters can be restricted to cover only the positive frequencies
%      and to be equivivalent to |dtwfbreal| by passing a `'real'` flag. 
%
%   `'freq'`(default),`'nat'`
%     The filters are ordered to produce subbands in the same order as 
%     |dtwfb| or |dtwfbreal| with the same flag.
%
%   Examples:
%   ---------
%
%   The following two examples create a multirate identity filterbank
%   using a duel-tree of depth 3:::
%
%      [g,a] = dtwfb2filterbank({'qshift3',3},'real');
%      filterbankfreqz(g,a,1024,'plot','linabs');
%
%   In the second example, the filterbank is identical to the full
%   wavelet tree:::
%
%      [g,a] = dtwfb2filterbank({'qshift3',3,'full'},'real');
%      filterbankfreqz(g,a,1024,'plot','linabs');
%
%   See also: dtwfbinit

complainif_notenoughargs(nargin,1,'DTWFB2FILTERBANK');

% Search for the 'real' flag
do_real = ~isempty(varargin(strcmp('real',varargin)));
if do_real
    %Remove the 'real' flag from varargin
    varargin(strcmp('real',varargin)) = [];
end

if ~isempty(varargin(strcmp('complex',varargin)))
    %Remove the 'complex' flag from varargin
    %It is not used elsewhere anyway
    varargin(strcmp('complex',varargin)) = [];
end


% Initialize the dual-tree
dtw = dtwfbinit({'strict',dualwt},varargin{:});

% Determine relation between the tree nodes
[wtPath, rangeLoc, rangeOut] = treeBFranges(dtw);
slice = ~cellfun(@isempty,rangeOut); % Limit to nodes with unconnected outputs
wtPath   = wtPath(slice); 
rangeLoc = rangeLoc(slice); 
rangeOut = rangeOut(slice);

% Multirate identity filters of the first tree
[g1,a] = nodesMultid(wtPath,rangeLoc,rangeOut,dtw);

% Multirate identity filters of the second tree
dtw.nodes = dtw.dualnodes;
g2 = nodesMultid(wtPath,rangeLoc,rangeOut,dtw);

if nargin>2
   % Return the filterbanks before doing the alignment
   info.g1 = cellfun(@(gEl) setfield(gEl,'h',gEl.h/2),g1,'UniformOutput',0);
   info.g2 = cellfun(@(gEl) setfield(gEl,'h',gEl.h/2),g2,'UniformOutput',0);
end

% Align filter offsets so they can be summed
for ii = 1:numel(g1)
   % Sanity checks
   assert(g1{ii}.offset<=0,sprintf('%s: Invalid wavelet filters.',upper(mfilename)));
   assert(g2{ii}.offset<=0,sprintf('%s: Invalid wavelet filters.',upper(mfilename)));
   
   % Insert zeros and update offsets
   offdiff = g1{ii}.offset-g2{ii}.offset;
   if offdiff>0
       g1{ii}.offset = g1{ii}.offset - offdiff;
       g1{ii}.h = [zeros(offdiff,1);g1{ii}.h(:)];
   elseif offdiff<0
       g2{ii}.offset = g2{ii}.offset + offdiff;
       g2{ii}.h = [zeros(-offdiff,1);g2{ii}.h(:)];
   end
   
   % Pad with zeros to a common length
   lendiff = numel(g1{ii}.h) - numel(g2{ii}.h);
   if lendiff~=0
       maxLen = max(numel(g1{ii}.h),numel(g2{ii}.h));
       g1{ii}.h = postpad(g1{ii}.h,maxLen);
       g2{ii}.h = postpad(g2{ii}.h,maxLen);
   end
end

% Filters covering the positive frequencies
g = cellfun(@(gEl,g2El) setfield(gEl,'h',(gEl.h+1i*g2El.h)/2),g1,g2,'UniformOutput',0);

% Mirror the filters when negative frequency filters are required too
if ~do_real
   gneg = cellfun(@(gEl,g2El) setfield(gEl,'h',(gEl.h-1i*g2El.h)/2),g1,g2,'UniformOutput',0);
   g = [g;gneg(end:-1:1)];
   a = [a;a(end:-1:1)];
end



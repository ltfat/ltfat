function [AF,BF]=wfbtbounds(wt,L)
%FWTBOUNDS Frame bounds of DWT
%   Usage: fcond=fwtbounds(wt,L);
%          [A,B]=fwtbounds(wt,L);
%
%   `wfbtbounds(wt,L)` calculates the ratio $B/A$ of the frame bounds
%   of the filterbank specified by *wt* for a system of length
%   *L*. The ratio is a measure of the stability of the system. 
%
%   `[A,B]=wfbtbounds(wt,L)` returns the lower and upper frame bounds
%   explicitly. 
%
%   See |wfbt| for explanation of parameter *wt*.
%
%   See also: wfbt, filterbankbounds


if nargin<2
  error('%s: Too few input parameters.',upper(mfilename));
end;

wt = wfbtinit({'strict',wt},'nat');

if L~=wfbtlength(L,wt)
    error(['%s: Specified length L is incompatible with the length of ' ...
           'the time shifts.'],upper(mfilename));
end;


for ii=1:numel(wt.nodes)
   a = wt.nodes{ii}.a;
   assert(all(a==a(1)),'%s: One of the basic wavelet filterbanks is not uniform.',upper(mfilename));
end

% Do the equivalent filterbank
[g,a] = wfbt2filterbank(wt);

maxGlen = max(cellfun(@(gEl) numel(gEl.h),g));
% maxa is also lcm(a)
maxa = max(a);

assert(L>=maxGlen,'%s: One of the filters is longer than L.',upper(mfilename));

gnew = {};

for m=1:numel(g)
   pk = maxa/a(m);
   
   for ii=0:pk-1
      gnew{end+1} = g{m};
      gnew{end+1}.offset = gnew{end+1}.offset-a(m)*ii;
   end
end

anew = maxa*ones(numel(gnew),1);


if nargout<2
   AF = filterbankbounds(gnew,anew,L);
elseif nargout == 2
   [AF, BF] = filterbankbounds(gnew,anew,L);
end
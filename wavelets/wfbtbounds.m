function [AF,BF]=wfbtbounds(wt,L)
%WFBTBOUNDS Frame bounds of WFBT
%   Usage: fcond=wfbtbounds(wt,L);
%          [A,B]=wfbtbounds(wt,L);
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

% Do the equivalent filterbank using multirate identity property
[gmultid,amultid] = wfbt2filterbank(wt);

% Check L
%maxGlen = max(cellfun(@(gEl) numel(gEl.h),gmultid));
%assert(L>=maxGlen,'%s: One of the filters is longer than L.',upper(mfilename));

% Do the equivalent uniform filterbank
[gu,au] = nonu2ufilterbank(gmultid,amultid);

if nargout<2
   AF = filterbankbounds(gu,au,L);
elseif nargout == 2
   [AF, BF] = filterbankbounds(gu,au,L);
end
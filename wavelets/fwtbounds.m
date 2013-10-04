function [AF,BF]=fwtbounds(w,J,L)
%FWTBOUNDS Frame bounds of DWT
%   Usage: fcond=fwtbounds(w,J,L);
%          [A,B]=fwtbounds(w,J,L);
%
%   `fwtbounds(w,a,L)` calculates the ratio $B/A$ of the frame bounds
%   of the filterbank specified by *w* and *J* for a system of length
%   *L*. The ratio is a measure of the stability of the system. 
%
%   `[A,B]=fwtbounds(w,J,L)` returns the lower and upper frame bounds
%   explicitly. 
%
%   See |fwt| for explanation of parameters *w* and *J*.
%
%   See also: fwt, filterbankbounds


if nargin<3
  error('%s: Too few input parameters.',upper(mfilename));
end;

if L~=fwtlength(L,w,J)
    error(['%s: Specified length L is incompatible with the length of ' ...
           'the time shifts.'],upper(mfilename));
end;

w = fwtinit({'strict',w});

assert(all(w.a==w.a(1)),'%s: Basic wavelet filterbank is not uniform.',upper(mfilename));

% There seem to be two possibilites:
% 1) Find the frame bounds of the equaivalent uniform filterbank. The
%    equivalent non-uniform filterbank can be created using fwt2filterbank,
%    non-uniform to uniform filterbank transform is described in the book:
%    Beyond Wavelets, chapter 10, Nonuniform filter banks: New results and
%    open problems by Akkarakaran and Vaidyanathan.
% 2) Use the recursive polyphase matrix multiplication algorithm from 
%    Stanhill, D.; Zeevi, Y. Y.: Frame analysis of wavelet-type filter banks 

%ad 1)
[g,a] = fwt2filterbank(w,J);

maxGlen = max(cellfun(@(gEl) numel(gEl.h),g));

assert(L>=maxGlen,'%s: One of the filters is longer than L.',upper(mfilename));

% amax is also lcm(a)
amax = max(a);

gnew = {};

for m=1:numel(g)
   pk = amax/a(m);
   
   for ii=0:pk-1
      gnew{end+1} = g{m};
      gnew{end+1}.offset = gnew{end+1}.offset-a(m)*ii;
   end
end

anew = amax*ones(numel(gnew),1);


if nargout<2
   AF = filterbankbounds(gnew,anew,L);
elseif nargout == 2
   [AF, BF] = filterbankbounds(gnew,anew,L);
end


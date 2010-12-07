function [AF,BF]=filterbankrealbounds(g,a,varargin);
%FILTERBANKBOUNDS  Frame bounds of filter bank for real signals only
%   Usage: fcond=filterbankbounds(g,a);
%          [A,B]=filterbankbounds(g,a);
%
%   FILTERBANKREALBOUNDS(g,a) calculates the ratio B/A of the frame bounds of
%   the filterbank specified by g and _a. The ratio is a measure of the
%   stability of the system.  Use this function on the common construction
%   where the filters in g only covers the positive frequencies.
%
%   [A,B]=FILTERBANKREALBOUNDS(g,a) returns the lower and upper frame bounds
%   explicitly.
%
%   See also: filterbank, filterbankdual
  
if nargin<2
  error('%s: Too few input parameters.',upper(mfilename));
end;

[M,longestfilter]=assert_filterbankinput(g,a);

[AF,BF]=filterbankbounds([g,cellfun(@conj,g,'UniformOutput',false)],a,varargin{:});

if nargout<2
  % Avoid the potential warning about division by zero.
  if AF==0
    AF=Inf;
  else
    AF=BF/AF;
  end;
end;
  

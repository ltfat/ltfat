function gd=filterbankrealdual(g,a,varargin);
%FILTERBANKDUAL  Dual filters of filterbank for real signals only 
%   Usage:  gd=filterbankdual(g,a);
%
%   FILTERABANKDUAL(g,a) computes the canonical dual filters of g for a
%   channel subsampling rate of _a (hop-size). The dual filters work only
%   for real-valued signals. Use this function on the common construction
%   where the filters in g only covers the positive frequencies.
%
%   The format of the filters g are described in the
%   help of FILTERANK.
%
%   To actually invert the output of a filterbank, use the dual filters
%   together with the IFILTERBANK function.
%
%   See also: filterbank, ifilterbank

if nargin<2
  error('%s: Too few input parameters.',upper(mfilename));
end;

[M,longestfilter]=assert_filterbankinput(g,a);

gd2=filterbankdual([g,cellfun(@conj,g,'UniformOutput',false)],a,varargin{:});

gd = gd2(M+1:2*M);


function L=filterbanklength(Ls,a);
%FILTERBANKLENGTH  Filterbank length from signal
%   Usage: L=filterbanklength(Ls,a);
%
%   `filterbanklength(Ls,a)` returns the length of a filterbank with
%   time shifts *a*, such that it is long enough to expand a signal of
%   length *Ls*.
%
%   If the filterbank length is longer than the signal length, the signal
%   will be zero-padded by |filterbank| or |ufilterbank|.
%
%   If instead a set of coefficients are given, call |filterbanklengthcoef|.
%
%   See also: filterbank, filterbanklengthcoef

if ~isnumeric(Ls)
  error('%s: Ls must be numeric.',upper(mfilename));
end;

if ~isscalar(Ls)
  error('%s: Ls must a scalar.',upper(mfilename));
end;

if ~isnumeric(a)
  error('%s: a must be numeric.',upper(mfilename));
end;

%if ~isvector(a)
%    
%end;

if any(a<=0)
      error('%s: "a" must consists of positive numbers only.',upper(mfilename));
end;

lcm_a=a(1);
for m=2:size(a,1)
  lcm_a=lcm(lcm_a,a(m,1));
end;

L=ceil(Ls/lcm_a)*lcm_a;

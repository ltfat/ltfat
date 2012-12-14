function L=nsdgtlength(Ls,a);
%NSDGTLENGTH  Nsdgt length from signal
%   Usage: L=nsdgtlength(Ls,a);
%
%   `nsdgtlength(Ls,a)` returns the length of an |nsdgt|_ with time shifts
%   *a*, such that it is long enough to expand a
%   signal of length *Ls*.
%
%   If the returned length is longer than the signal length, the signal
%   will be zero-padded by |nsdgt|_ or |unsdgt|_.
%
%   If instead a set of coefficients are given, call |nsdgtlengthcoef|_.
%
%   See also: nsdgt, nsdgtlengthcoef

if ~isnumeric(Ls)
  error('%s: Ls must be numeric.',upper(mfilename));
end;

if ~isscalar(Ls)
  error('%s: Ls must a scalar.',upper(mfilename));
end;

if ~isnumeric(a)
  error('%s: a must be numeric.',upper(mfilename));
end;

if ~isvector(a) || any(a<0)
  error('%s: "a" must be a vector of non-negative numbers.',upper(mfilename));
end;

L=sum(a);

if Ls>L
    error('%s: The signal must have fewer than %i samples.',upper(mfilename),L);
end;
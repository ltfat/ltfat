function L=nsdgtlength(Ls,a);
%NSDGTLENGTH  NSDGT length from signal
%   Usage: L=nsdgtlength(Ls,a);
%
%   `nsdgtlength(Ls,a)` returns the length of an |nsdgt| with time shifts
%   *a*, such that it is long enough to expand a
%   signal of length *Ls*.
%
%   If the returned length is longer than the signal length, the signal
%   will be zero-padded by |nsdgt| or |unsdgt|.
%
%   If instead a set of coefficients are given, call |nsdgtlengthcoef|.
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
    error('%s: The signal must have at most %i samples.',upper(mfilename),L);
end;

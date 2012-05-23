function L=nonsepdgtlengthsignal(Ls,a,M,s);
%NONSEPDGTLENGTHSIGNAL  Non-separable DGT length from signal
%   Usage: L=nonsepdgtlengthsignal(Ls,a,M,s);
%
%   `nonsepdgtlengthsignal(Ls,a,M,s)` returns the length of a Gabor system
%   on a non-separable lattice that is long enough to expand a signal of
%   length *Ls*. Please see the help on |nonsepdgt|_ for an explanation
%   of the parameters *a*, *M* and *s*.
%
%   If the returned length is longer than the signal length, the signal
%   will be zero-padded by |nonsepdgt|_.
%
%   See also: nonsepdgt

if ~isscalar(M)
  error('%s: M must be a scalar',upper(mfilename));
end;

if ~isscalar(a)
  error('%s: a must be a scalar',upper(mfilename));
end;

if rem(M,1)~=0
  error('%s: M must be an integer',upper(mfilename));
end;

if rem(a,1)~=0
  error('%s: a must be an integer',upper(mfilename));
end;

if ~isnumeric(s) || ~isvector(s) || length(s)~=2
    error('%s: s must be a vector of length 2.',upper(mfilename));
end;

if ~isnumeric(Ls)
  error('%s: Ls must be numeric.',upper(mfilename));
end;

if ~isscalar(Ls)
  error('%s: Ls must a scalar.',upper(mfilename));
end;

Lsmallest=lcm(a*s(2),M*s(2));

L=ceil(Ls/Lsmallest)*Lsmallest;



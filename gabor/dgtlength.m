function [L,tfr]=dgtlength(Ls,a,M,lt);
%DGTLENGTH  DGT length from signal
%   Usage: L=dgtlength(Ls,a,M,s);
%
%   `dgtlength(Ls,a,M)` returns the length of a Gabor system that is long
%   enough to expand a signal of length *Ls*. Please see the help on
%   |dgt|_ for an explanation of the parameters *a* and *M*.
%
%   If the returned length is longer than the signal length, the signal
%   will be zero-padded by |dgt|_.
%
%   `dgtlength(Ls,a,M,lt)` does the same for a non-separable lattice.
%
%   See also: dgt

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

if ~isnumeric(Ls)
    error('%s: Ls must be numeric.',upper(mfilename));
end;

if ~isscalar(Ls)
    error('%s: Ls must a scalar.',upper(mfilename));
end;

if nargin<4

    Lsmallest=lcm(a,M);    

else

    if ~isnumeric(lt) || ~isvector(lt) || length(lt)~=2
        error('%s: lt must be a vector of length 2.',upper(mfilename));
    end;
    

    if (mod(lt(2),1)>0) || lt(2)<=0
        error('%s: lt(2) must be a positive integer.',upper(mfilename));
    end;
    
    if (mod(lt(1),1)>0) || lt(1)<0 || lt(1)>=lt(2)
        error(['%s: lt(1)=%i must be a positive integer that is larger than 0 but ' ...
               'smaller than lt(2)=%i.'],upper(mfilename),lt(1),lt(2));
    end;

    if lt(1)==0 && lt(2)~=1
        error('%s: The rectangular lattice can only be specified by lt=[0 1].',upper(mfilename));
    end;

    if gcd(lt(1),lt(2))>1
        error('%s: lt(1)/lt(2) must be an irriducible fraction.',upper(mfilename));
    end;

    Lsmallest=lcm(a,M)*lt(2);

end;

L=ceil(Ls/Lsmallest)*Lsmallest;



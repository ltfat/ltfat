function [L,tfr]=dwiltlength(Ls,M);
%DWILTLENGTH  DWILT/WMDCT length from signal
%   Usage: L=dwiltlength(Ls,M);
%
%   `dwiltlength(Ls,M)` returns the length of a Wilson / WMDCT system with
%   *M* channels system is long enough to expand a signal of length
%   *Ls*. Please see the help on |dwilt| or |wmdct| for an explanation of the
%   parameter *M*.
%
%   If the returned length is longer than the signal length, the signal will
%   be zero-padded by |dwilt| or |wmdct|.
%
%   A valid transform length must be divisable by *2M*. This
%   means that the minumal admissable transform length is ::
%
%     Lsmallest = 2*M;
%
%   and all valid transform lengths are multipla of *Lsmallest*
%
%   See also: dwilt, wmdct

complainif_argnonotinrange(nargin,2,2,mfilename);

if ~isnumeric(M) || ~isscalar(M)
  error('%s: M must be a scalar',upper(mfilename));
end;

if rem(M,1)~=0 || M<=0
  error('%s: M must be a positive integer',upper(mfilename));
end;

if ~isnumeric(Ls)
    error('%s: Ls must be numeric.',upper(mfilename));
end;

if ~isscalar(Ls)
    error('%s: Ls must a scalar.',upper(mfilename));
end;

Lsmallest=2*M;

L=ceil(Ls/Lsmallest)*Lsmallest;

b=L/(2*M);
tfr=M/b;



function f=prect(L,n)
%PRECT   Periodic rectangle
%   Usage:  f=prect(L,n);
%
%   `psinc(L,n)` computes the periodic rectangle (or square) function of
%   length *L* supported on *n* samples. The |dft| of the periodic
%   rectangle function in the periodic sinc function, |psinc|.
%
%   * If *n* is odd, the output will be supported on exactly *n* samples
%     centered around the first sample.
%
%   * If *n* is even, the output will be supported on exactly *n+1* samples
%     centered around the first sample. The function value on the two
%     samples on the edge of the function will have half the magnitude of
%     the other samples.
%
%   Examples:
%   ---------
%
%   This figure displays an odd length periodic rectangle:::
%
%     stem(prect(30,11));
%     ylim([-.2 1.2]);
%
%   This figure displays an even length periodic rectangle. Notice the
%   border points:::
%
%     stem(prect(30,12));
%     ylim([-.2 1.2]);
%
%   See also: psinc

complainif_argnonotinrange(nargin,2,2,mfilename);

if ~(numel(L)==1) || ~(isnumeric(L)) || mod(L,1)~=0 || L<=0
    error('%s: L has to be a positive integer.',upper(mfilename));
end;

if ~(numel(n)==1) || ~(isnumeric(L)) || mod(n,1)~=0 || n<=0
    error('%s: n has to be a positive integer.',upper(mfilename));
end;

f=pbspline(L,0,n);


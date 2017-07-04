function f=psinc(L,n)
%PSINC   Periodic Sinc function (Dirichlet function)
%   Usage:  f=psinc(L,n);
%
%   `psinc(L,n)` computes the periodic Sinc function of length *L* with
%   $n-1$ local extrema. The |dft| of the periodic Sinc function is the
%   periodic rectangle, |prect|, of length *n*.
%
%   Examples:
%   ---------
%
%   This figure displays a the periodic sinc function with 6 local extremas:::
%
%     plot(psinc(30,7));
%
%   See also: prect

complainif_argnonotinrange(nargin,2,2,mfilename);

if ~(numel(L)==1) || ~(isnumeric(L)) || mod(L,1)~=0 || L<=0
    error('%s: L has to be a positive integer.',upper(mfilename));
end;

if ~(numel(n)==1) || ~(isnumeric(L)) || mod(n,1)~=0 || n<=0
    error('%s: n has to be a positive integer.',upper(mfilename));
end;

x=(2*pi*(0:L-1)/L).';

n_odd = n-(1-mod(n,2));

f = sin(n_odd.*x./2)./(n_odd.*sin(x./2));

f(1)  = 1;

if (mod(n,2))==0;
    f = f+cos(x*n/2)/n_odd;
end;


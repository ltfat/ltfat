function d=nsgabframediag(g,a,M)
%NSGABFRAMEDIAG  Diagonal of Gabor frame operator
%   Usage:  d=nsgabframediag(g,a,M);
%
%   Input parameters:
%         g     : Window function.
%         a     : Length of time shift.
%         M     : Number of channels.
%   Output parameters:
%         d     : Diagonal stored as a column vector
%
%   `nsgabframediag(g,a,M)` computes the diagonal of the non-stationary
%   Gabor frame operator with respect to the window *g* and parameters *a*
%   and *M*. The diagonal is stored as a column vector of length `L=sum(a)`.
%
%   The diagonal of the frame operator can for instance be used as a
%   preconditioner.
%
%   See also: nsdgt

if nargin<3
  error('%s: Too few input parameters.',upper(mfilename));
end;

L=sum(a);

timepos=cumsum(a)-a(1);

N=length(a);

[g,info]=nsgabwin(g,a,M);

a=info.a;
M=info.M;

d=zeros(L,1,assert_classname(g{1}));
for ii=1:N
    shift=floor(length(g{ii})/2);
    temp=abs(circshift(g{ii},shift)).^2*M(ii);
    tempind=mod((1:length(g{ii}))+timepos(ii)-shift-1,L)+1;
    d(tempind)=d(tempind)+temp;
end



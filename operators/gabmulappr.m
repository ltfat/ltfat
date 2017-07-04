function [sym,lowb,upb]=gabmulappr(T,p2,p3,p4,p5);
%GABMULAPPR  Best Approximation by a Gabor multiplier
%   Usage:  sym=gabmulappr(T,a,M);
%           sym=gabmulappr(T,g,a,M);
%           sym=gabmulappr(T,ga,gs,a,M);
%           [sym,lowb,upb]=gabmulappr( ... );
%
%   Input parameters:
%         T     : matrix to be approximated
%         g     : analysis/synthesis window
%         ga    : analysis window
%         gs    : synthesis window
%         a     : Length of time shift.
%         M     : Number of channels.
%
%   Output parameters:
%         sym   : symbol
%
%   `sym=gabmulappr(T,g,a,M)` calculates the best approximation of the given
%   matrix *T* in the Frobenius norm by a Gabor multiplier determined by the
%   symbol *sym* over the rectangular time-frequency lattice determined by
%   *a* and *M*.  The window *g* will be used for both analysis and
%   synthesis.
%
%   `gabmulappr(T,a,M)` does the same using an optimally concentrated, tight
%   Gaussian as window function.
%
%   `gabmulappr(T,gs,ga,a)` does the same using the window *ga* for analysis
%   and *gs* for synthesis.
%
%   `[sym,lowb,upb]=gabmulappr(...)` additionally returns the lower and
%   upper Riesz bounds of the rank one operators, the projections resulting
%   from the tensor products of the analysis and synthesis frames.
%
%   See also: framemulappr
%
%   Demos: demo_gabmulappr
%
%   References: doetor09 ba07 ba06 fehakr04

% AUTHOR    : Monika Döerfler
% REFERENCE : REF_GABMULAPPR
% TESTING   : TEST_GABMULAPPR
  
complainif_argnonotinrange(nargin,3,5,mfilename);

L=size(T,1);

if size(T,2)~=L
  error('T must be square.');
end;

if nargin==3
  % Usage: sym=gabmulappr(T,a,M);
  a=p2;
  M=p3;
  ga=gabtight(a,M,L);
  gs=ga;
end;

if nargin==4
  % Usage: sym=gabmulappr(T,g,a,M);
  ga=p2;
  gs=p2;
  a=p3;
  M=p4;
end;
  
if nargin==5
  % Usage: sym=gabmulappr(T,ga,gm,a,M);
  ga=p2;
  gs=p3;
  a=p4;
  M=p5;
end;

if size(ga,2)>1
  if size(ga,1)>1
    error('Input g/ga must be a vector');
  else
    % ga was a row vector.
    ga=ga(:);
  end;
end;

if size(gs,2)>1
  if size(gs,1)>1
    error('Input g/gs must be a vector');
  else
    % gs was a row vector.
    gs=gs(:);
  end;
end;

N=L/a;
b=L/M;

Vg=dgt(gs,ga,1,L);

s=spreadfun(T);

A=zeros(N,M);
V=zeros(N,M);
for k=0:b-1 
  for l=0:a-1
    A = A+ s(l*N+1:(l+1)*N,k*M+1:k*M+M).*conj(Vg(l*N+1:(l+1)*N,k*M+1:k*M+M));
    V = V+abs(Vg(l*N+1:(l+1)*N,k*M+1:k*M+M)).^2;
  end;
end;

if nargout>1
  lowb = min(V(:));
  upb  = max(V(:));
end;

SF1=A./V;

SF=zeros(N,M);
jjmod=mod(-M:-1,M)+1;
iimod=mod(-N:-1,N)+1;
SF=SF1(iimod,jjmod);

sym=b*dsft(SF)*sqrt(M)/sqrt(N);

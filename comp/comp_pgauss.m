function [g]=comp_pgauss(L,w,center)
%COMP_PGAUSS  Sampled, periodized Gaussian.
%   
%   Computational routine: See help on PGAUSS.
%
%   center=0  gives whole-point centered function.
%   center=.5 gives half-point centered function.
%
%   Does not check input parameters, do not call this
%   function directly.

%   AUTHOR : Peter Soendergaard.

%   AUTHOR : Peter Soendergaard.
%   TESTING: OK
%   REFERENCE: OK

if prod(size(center))==1
  c_t=center;
  c_f=0;
else
  c_t=center(1);
  c_f=center(2);
end;

% c_t - time centering
% c_f - frequency centering

g=zeros(L,1);

if L==0
  return;
end;

sqrtl=sqrt(L);
safe=4;

% Outside the interval [-safe,safe] then exp(-pi*x.^2) is numerically zero.
nk=ceil(safe/sqrt(L/sqrt(w)));

lr=(0:L-1).';
for k=-nk:nk  
  % The following line includes a frequency center. This is not exposed in the
  % current code 
  %g=g+exp(-pi*((lr+c_t)/sqrtl-k*sqrtl).^2/w-2*pi*i*c_f*(lr/L-k));

  g=g+exp(-pi*((lr+c_t)/sqrtl-k*sqrtl).^2/w);
end;

% Normalize it exactly.
g=g/norm(g);

% This normalization is only approximate
%g=g*(w*L/2)^(-.25);



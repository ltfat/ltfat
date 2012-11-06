function x=modcent(x,r);
%MODCENT  Centered modulo
%   Usage:  y=modcent(x,r);
%
%   `modcent(x,r)` computes the modulo of *x* in the range $[-r/2,r/2[$.
%
%   As an example, to compute the modulo of *x* in the range $[-\pi,\pi[$ use
%   the call::
%
%     y = modcent(x,2*pi);

x=mod(x,r);  
idx=x>r/2;
x(idx)=x(idx)-r;

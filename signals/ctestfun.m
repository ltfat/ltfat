function [ftest]=ctestfun(L)
%CTESTFUN  Complex 1-D test function
%   Usage:  ftest=ctestfun(L);
%
%   `ctestfun(L)` returns a test signal consisting of a superposition of a
%   chirp and an indicator function.

ftest=zeros(L,1);

sp=round(L/4);
lchirp=round(L/2);
ftest(sp+1:sp+lchirp)=exp(2*i*linspace(0,2*pi*sqrt(lchirp)/10,lchirp).^2)';

s=round(L*7/16);
l=round(L/16);
ftest(s:s+l)=ftest(s:s+l)+ones(l+1,1);


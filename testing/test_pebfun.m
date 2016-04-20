clear all
close all
clc

w = [-3,-1,1,3];
a = 25;
M = 31;
inc = 10;
L = 1e6; L = dgtlength(L,a,M);
width = M;

g = pebfun(w,L,width);
tic
gamma = pebfundual(w,a,M,L,width,inc);
toc

tic
gammaLTFAT = gabdual(g,a,M,L);
toc

f = randn(1,L);

display('gabdualnorm "new approach":')
[scal,err] = gabdualnorm(g,gamma,a,M,L)
display('gabdualnorm LTFAT:')
[scal,err] = gabdualnorm(g,gammaLTFAT,a,M,L)

fr = idgt(dgt(f,g,a,M),gamma,a,numel(f));
fr = fr(:)';
display('reconstruction error "new approach":')
norm(f-fr)/norm(f)

fr = idgt(dgt(f,g,a,M),gammaLTFAT,a,numel(f));
fr = fr(:)';
display('reconstruction error LTFAT:')
norm(f-fr)/norm(f)

% norm(gamma-gammaLTFAT)
% plot([1:L],[fftshift(gamma)';fftshift(gammaLTFAT)'])
% plot([1:L],(fftshift(gamma)'-fftshift(gammaLTFAT)'))
% axis([L/2-1e3,L/2+1e3,ylim])
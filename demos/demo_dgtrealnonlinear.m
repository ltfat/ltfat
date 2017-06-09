function demo_dgtrealnonlinear
f = randn(2*512,1);
a = 16;M=512;
c = dgt(f,'hann',a,M);
p = 0;
val = exp(i*pi/3);
pos = 21;

c(:) = 0;
c(1,1) = 1;   
kern =  phasekernfi(dgt(idgt(c,'hann',a),'hann',a,M),pos-1,a,M);

K = sqrt(2/(1 + real(exp(i*2*angle(val))*conj(kern(3+2*p,1)))) )


c(:) = 0;
c([2+p],pos) = val;   
%c([end-p],pos) = conj(val);
fhat1 = idgt(c,'hann',a);
%figure(1);plot([real(fhat1),imag(fhat1)]);shg

c(:) = 0;
c([2+p],pos) = val; 
c([end-p],pos) = conj(val);
fhat2 = idgt(c,'hann',a);
figure(1);plot([real(fhat1),real(fhat2)]);shg

Ktrue = 1/norm(real(fhat2/2))
%norm(fhat2)^2


% c2 = dgt(fhat1,'hann',a,M);
% 
% figure(1); plotdgt(c,a,'linabs','tc');
% figure(2); plotdgt(imag(c2),a,'lin','tc');
% figure(3); plotdgt(real(c2),a,'lin','tc');
% 
% c3 = val*circshift(kern,[p+1,pos-1]) + conj(val)*circshift(kern,[-(p+1),pos-1]);
% 
% figure(4); plotdgt(imag(c3),a,'lin','tc');
% figure(5); plotdgt(real(c3),a,'lin','tc');

function kernm = phasekernfi(kern,n,a,M)
    kernh = size(kern,1);
    l = -2*pi*n*fftindex(kernh)*a/M;
    kernm = bsxfun(@times,kern,exp(1i*l));
f = gspi;

J = 10;
w = fwtinit({'db',10});
[h,a] = multid(w,J);
H = freqzfb(h,filterbanklength(length(f),a));
c = wtfft(f,H,a);

figure(1);
plotfwt(c);
figure(2);
plot(f);
%figure(3);
%freqzfb(h,nextpow2(length(f)));

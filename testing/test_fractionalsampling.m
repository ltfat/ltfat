f=[1;zeros(299,1)];
g=blfilter('hanning',.2,.2);

oo=pfilt(f,g,[3 2]);

subplot(2,1,1);
magresp(g,'dynrange',90);

subplot(2,1,2);
plotfft(fft(oo),'dynrange',90);

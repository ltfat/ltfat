%f = randn(10000,1);
f = zeros(100,1);
f(1) = 1;

wt =  {'qshift3',2,'first','ana:symorth3'};

[g,a,info] = dtwfbreal2filterbank( wt,'nat');
c = filterbank(f,g,a);

wt2 =  dtwfbinit(wt);

figure(1);
clf;
F =filterbankfreqz(g,a,2048,'linabs','plot');
subplot(2,1,1);
plot(abs(fftshift(F)));
subplot(2,1,2);
plot(unwrap(angle(fftshift(F))));

cd = dtwfbreal(f,wt);


figure(2);
clf;
C = linspecer(2);


L = max(cellfun(@(gEl) numel(gEl.h),g))*2;
%n = -floor(L/2)+1:floor(L/2);
n = 0:L-1;
for ii=1:numel(g)
    clf;
    H1imag = (fftshift(fft(circshift(postpad(info.g1{ii}.h,L),info.g1{ii}.offset))));
    h1imag = (circshift(postpad(info.g1{ii}.h,L),info.g1{ii}.offset));
    h = stem(n,imag(H1imag));
    hold on;
    set(h,'Color',C(1,:));
    drawnow;
    
    H2imag = (fftshift(fft(circshift(postpad(info.g2{ii}.h,L),info.g2{ii}.offset))));
    h2imag = (circshift(postpad(info.g2{ii}.h,L),info.g2{ii}.offset));
    h = stem(n,imag(H2imag));
    set(h,'Color',C(2,:));
    drawnow;
    
    stem(n,abs(H1imag+i*H2imag),'k')
    
    pause(0.5);
end
hold off;

nodId = 1;
g1=wt2.nodes{nodId}.g;
g2=wt2.dualnodes{nodId}.g;

for ii=1:numel(g1)
    figure(8);clf;
    h1 = (circshift(postpad(g1{ii}.h,L),g1{ii}.offset));
    H1imag = (angle(fft(circshift(postpad(g1{ii}.h,L),g1{ii}.offset))));
     
    h = stem(n/L*2*pi,H1imag);
    hold on;
    set(h,'Color',C(1,:));
    drawnow;
    
    h2 = (circshift(postpad(g2{ii}.h,L),g2{ii}.offset));
    H2imag = (angle(fft(circshift(postpad(g2{ii}.h,L),g2{ii}.offset))));
     
    h = stem(n/L*2*pi,H2imag);
    hold on;
    set(h,'Color',C(2,:));
    drawnow;
    
end
    


figure(4);
stem(imag([(cell2mat(c)),cell2mat(cd)]))
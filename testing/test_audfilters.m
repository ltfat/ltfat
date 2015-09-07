% Testing script for audfilters

[f,fs]=greasy;  % Get the test signal
L=length(f);


[g,a,fc]=audfilters(fs,L,0,fs/2,'mel','fractional');
c=filterbank(f,{'realdual',g},a);
r=2*real(ifilterbank(c,g,a));
norm(f-r)

% Plot the response
figure(1);
subplot(2,1,1);
R=filterbankresponse(g,a,L,fs,'real','plot');

subplot(2,1,2);
semiaudplot(linspace(0,fs/2,L/2+1),R(1:L/2+1));
ylabel('Magnitude');

% Plot frequency responses of individual filters
gd=filterbankrealdual(g,a,L);
figure(2);
subplot(2,1,1);
filterbankfreqz(gd,a,L,fs,'plot','linabs','posfreq');

subplot(2,1,2);
filterbankfreqz(g,a,L,fs,'plot','linabs','posfreq');
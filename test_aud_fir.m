clear all;

len = 5;
fs = 16000;
Ls = fs*len;

fmin = 100;
fmax = fs/2;
num_filters = 128;
mul = 32;

[g, a, fc, L] = audfilters(fs,Ls,'M',num_filters,'uniform','mel','bwmul',mul,'redmul',1);

[g_fir, a_fir, fc_fir, L_fir] = audfilters_fir2(fs,Ls,'M',num_filters,'uniform','mel','bwmul',mul,'redmul',1);

%implicitly calculated in audfilters
%L = filterbanklength(Ls,a);
%L_fir = filterbanklength(Ls,a_fir);

figure
gf = filterbankresponse(g,a,L,'plot','real');
%figure
hold on
gf_fir = filterbankresponse(g_fir,a_fir,L_fir,'plot','real');

f_test = rand(Ls,1);

c=ufilterbank(f_test,g_fir,a_fir);
r=2*real(ifilterbank(c,{'realdual',g_fir},a_fir));

sprintf('Reconstruction error FIR FB: %d', norm(f_test-r))

fb = filterbankrealbounds(g_fir,a_fir,L);

sprintf('Framebound ratio FIR FB: %d', fb)

figure;
gf = filterbankfreqz(g_fir,a_fir,L,'plot');
%the filters "only" go down to -300 dB (probably due to truncation)
ylim([-300 100])

%plotfilterbank(c,a,fc,fs,90,'audtick');

gt = filterbanktight(g_fir,a_fir,L);
c=filterbank(f_test,gt,a_fir);
r=2*real(ifilterbank(c,{'realdual',gt},a_fir));

sprintf('Reconstruction error FIR FB tight: %d', norm(f_test-r))

fbtight = filterbankrealbounds(gt,a_fir,L);

sprintf('Framebound ratio FIR FB tight: %d', fbtight)

figure;
gf = filterbankfreqz(gt,a_fir,L,'plot');
ylim([-300 100])
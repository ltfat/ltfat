len = 5;
fs = 16000;
Ls = fs*len;

fmin = 100;
fmax = fs/2;
num_filters = 128;
mul = 32;

[g, a, fc, L] = audfilters(fs,Ls,'M',num_filters,'uniform','mel','bwmul',mul,'redmul',1);

[g_fir, a_fir, fc_fir, L] = audfilters_fir2(fs,Ls,'M',num_filters,'uniform','mel','bwmul',mul,'redmul',1);

L = filterbanklength(Ls,a);
L_fir = filterbanklength(Ls,a_fir);

figure
gf = filterbankresponse(g,a,L,'plot','real');
figure
gf_fir = filterbankresponse(g_fir,a_fir,L_fir,'plot','real');

%c=filterbank(rand(Ls,1),g,a);

%plotfilterbank(c,a,fc,fs,90,'audtick');

%gt = filterbanktight(g,a,L);

%filterbankrealbounds(gt,a,L)
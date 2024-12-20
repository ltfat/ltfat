len = 5;
fs = 16000;
Ls = 16000*5;

fmin = 100;
fmax = fs/2;
num_filters = 512;
red = 8;

[g, a, fc, L] = audfilters_fir(fs,Ls,'M',num_filters,'uniform','mel','bwmul',16,'redtar',3.98);

%gt = filterbanktight(g,a,L);
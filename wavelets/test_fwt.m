function test_failed=test_fwt
%%
f = randn(17897,1);
J = 1;
% Examples


w = waveletfb('algmband',1);
[c] = fwt(f,w,J);
fhat = ifwt(c,w,J,length(f));
stem([f,fhat]);
%%
w = waveletfb('db',8);
[c] = fwt(f,w,J);
fhat = ifwt(c,w,J,length(f));
stem([f,fhat]);
%%

[h,g] = wfilt_db(8);
[c] = fwt(f,h,J);
fhat = ifwt(c,g,J,length(f));
stem([f,fhat]);
%%

w = waveletfb('db',8,'undec');
[c] = fwt(f,w,J);
fhat = ifwt(c,w,J,length(f));
stem([f,fhat]);

w = waveletfb('db',8,'undec','sym');
[c] = fwt(f,w,J);
fhat = ifwt(c,w,J,length(f));
stem([f,fhat]);


w = waveletfb('db',8,'undec','sym');
[c] = fwt(f,w,J,'per');
fhat = ifwt(c,w,J,length(f),'per');
stem([f,fhat]);

w = waveletfb('db',8,'undec','sym');
[c] = fwt(f,w,J);
fhat = ifwt(c,w,J);
stem([f,fhat]);


w = waveletfb('db',8);
[c] = fwt(f,w,J,'undec','sym');
fhat = ifwt(c,w,J,'undec','sym');
stem([f,fhat]);


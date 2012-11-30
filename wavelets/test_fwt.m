function test_failed=test_fwt

f = randn(17897,1);
J = 10;
% Examples

w = waveletfb('db',8);
[c] = fwt(f,w,J);
fhat = ifwt(c,w,length(f));
stem([f,fhat]);

[h,g] = dbfilt(8);
[c] = fwt(f,h,J);
fhat = ifwt(c,g,length(f));
stem([f,fhat]);

[h,g] = dbfilt(8);
[c,Lc] = fwt(f,h,J);
fhat = ifwt(c,Lc,g,length(f));
stem([f,fhat]);


w = waveletfb('db',8,'undec');
[c] = fwt(f,w,J);
fhat = ifwt(c,w,length(f));
stem([f,fhat]);

w = waveletfb('db',8,'undec','sym');
[c] = fwt(f,w,J);
fhat = ifwt(c,w,length(f));
stem([f,fhat]);


w = waveletfb('db',8,'undec','sym');
[c] = fwt(f,w,J,'per');
fhat = ifwt(c,w,length(f),'per');
stem([f,fhat]);

w = waveletfb('db',8,'undec','sym');
[c] = fwt(f,w,J);
fhat = ifwt(c,w);
stem([f,fhat]);


w = waveletfb('db',8);
[c] = fwt(f,w,J,'undec','sym');
fhat = ifwt(c,w,'undec','sym');
stem([f,fhat]);


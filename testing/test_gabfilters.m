function test_failed = test_gabfilters
test_failed = 0;

M = 512; a = 128; L = 10*512; g = randn(M,1);
N = L/a;
M2 = floor(M/2) + 1;

f = randn(L,1);

c1 = dgt(f,g,a,M,'timeinv');
[gfb,afb] = gabfilters(L,g,a,M,'complex');
c2 = ufilterbank(f,gfb,afb);
res = norm(c1 - c2.');

[test_failed,fail]=ltfatdiditfail(res,test_failed);
fprintf('GABFILTERS DGT     UFILTERBANK TIME %s\n',fail);

c3 = filterbank(f,gfb,afb);
c3 = reshape(cell2mat(c3),N,M);
res = norm(c1 - c3.');

[test_failed,fail]=ltfatdiditfail(res,test_failed);
fprintf('GABFILTERS DGT      FILTERBANK TIME %s\n',fail);

gfft = fir2long(g,L);
g = ifft(fir2long( gfft , L));

c1 = dgt(f,g,a,M,'timeinv');
[gfb,afb] = gabfilters(L,gfft,a,M,'freq','complex');

c2 = ufilterbank(f,gfb,afb);
res = norm(c1 - c2.');

[test_failed,fail]=ltfatdiditfail(res,test_failed);
fprintf('GABFILTERS DGT     UFILTERBANK FREQ %s\n',fail);

c3 = filterbank(f,gfb,afb);
c3 = reshape(cell2mat(c3),N,M);
res = norm(c1 - c3.');

[test_failed,fail]=ltfatdiditfail(res,test_failed);
fprintf('GABFILTERS DGT      FILTERBANK FREQ %s\n',fail);

g = randn(M,1);
c1 = dgtreal(f,g,a,M,'timeinv');
[gfb,afb] = gabfilters(L,g,a,M);
c2 = ufilterbank(f,gfb,afb);
res = norm(c1 - c2.');

[test_failed,fail]=ltfatdiditfail(res,test_failed);
fprintf('GABFILTERS DGTREAL UFILTERBANK TIME %s\n',fail);

c3 = filterbank(f,gfb,afb);
c3 = reshape(cell2mat(c3),N,M2);
res = norm(c1 - c3.');

[test_failed,fail]=ltfatdiditfail(res,test_failed);
fprintf('GABFILTERS DGTREAL  FILTERBANK TIME %s\n',fail);


gfftlong = fir2long(g,L);
gfft = gfftlong + involute(gfftlong);
g = real(ifft(fir2long( gfft , L)));

c1 = dgtreal(f,g,a,M,'timeinv');
[gfb,afb] = gabfilters(L,gfft,a,M,'freq');
c2 = ufilterbank(f,gfb,afb);

res = norm(c1 - c2.');

[test_failed,fail]=ltfatdiditfail(res,test_failed);
fprintf('GABFILTERS DGTREAL UFILTERBANK FREQ %s\n',fail);

c3 = filterbank(f,gfb,afb);
c3 = reshape(cell2mat(c3),N,M2);
res = norm(c1 - c3.');

[test_failed,fail]=ltfatdiditfail(res,test_failed);
fprintf('GABFILTERS DGTREAL  FILTERBANK FREQ %s\n',fail);
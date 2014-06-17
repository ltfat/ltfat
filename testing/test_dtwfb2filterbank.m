function test_failed = test_dtwfb2filterbank
test_failed = 0;

L = 1000;

dualwt = {'ana:qshift1',3};

[g,a,info] = dtwfbreal2filterbank( dualwt);

G = filterbankfreqz(g,a,L);

Greal = filterbankfreqz(info.g1,a,L);
Ghilb = filterbankfreqz(info.g2,a,L);

G2 = Greal+1i*Ghilb;


res = norm(G(:)-G2(:));

[test_failed,fail]=ltfatdiditfail(res,test_failed);
s=sprintf(['DUAL-TREE  %0.5g %s'],res,fail);
disp(s)

figure(1);
plot([abs(G)]);

figure(2);
plot([abs(G2)]);

function test_failed=test_rangecompress
%TEST_RANGECOMPRESS Test range compression and expansion

test_failed=0;

disp(' ===============  TEST_RANGECOMPRESS ================');

x=tester_crand(5,11);

y=rangecompress(x,'mulaw');
x_r=rangeexpand(y,'mulaw');

res=norm(x-x_r,'fro');

[test_failed,fail]=ltfatdiditfail(res,test_failed);
fprintf(['RANGECOMPRESS MULAW %0.5g %s\n'],res,fail);

y=rangecompress(x,'alaw');
x_r=rangeexpand(y,'alaw');

res=norm(x-x_r,'fro');

[test_failed,fail]=ltfatdiditfail(res,test_failed);
fprintf(['RANGECOMPRESS  ALAW %0.5g %s\n'],res,fail);



function test_failed=test_gga
disp(' ===============  TEST_GGA ================');
test_failed=0;
% Goertzel algorithm has a bad precision for larger L.
global LTFAT_TEST_TYPE;
tolerance = 2e-10;
if strcmpi(LTFAT_TEST_TYPE,'single')
   tolerance = 3e-4;
end
L = 36;
W = 17;
f=tester_crand(L,W); 



res = norm(fft(cast(f,'double'))-gga(f,0:L-1));
[test_failed,fail]=ltfatdiditfail(res,test_failed,tolerance);
fprintf('RES 1 cols: L:%3i, W:%3i %s\n',L,W,fail);

res = norm(fft(cast(f,'double'),[],2)-gga(f,0:W-1,2));
[test_failed,fail]=ltfatdiditfail(res,test_failed,tolerance);
fprintf('RES 1 rows: L:%3i, W:%3i %s\n',L,W,fail);


 res = norm(fft(cast(f,'double'),2*L)-gga(f,0:0.5:L-0.5));
 [test_failed,fail]=ltfatdiditfail(res,test_failed,tolerance);
 fprintf('RES 1/2 cols: L:%3i, W:%3i %s\n',L,W,fail);
 
 res = norm(fft(cast(f,'double'),5*L)-gga(f,0:1/5:L-1/5));
 [test_failed,fail]=ltfatdiditfail(res,test_failed,tolerance);
 fprintf('RES 1/5 cols: L:%3i, W:%3i %s\n',L,W,fail);
function test_failed=test_gga
disp(' ===============  TEST_GGA ================');
test_failed=0;
% Goertzel algorithm has a bad precision for larger L.
global LTFAT_TEST_TYPE;
tolerance = 2e-10;
if strcmpi(LTFAT_TEST_TYPE,'single')
   tolerance = 9e-3;
end
L = 36;
W = 17;
f=tester_crand(L,W); 



res = norm(fft(cast(f,'double'))-gga(f,linspace(0,1-1/L,L)));
[test_failed,fail]=ltfatdiditfail(res,test_failed,tolerance);
fprintf('RES 1 cols: L:%3i, W:%3i %s\n',L,W,fail);

res = norm(fft(cast(f,'double'),[],2)-gga(f,linspace(0,1-1/W,W),[],2));
[test_failed,fail]=ltfatdiditfail(res,test_failed,tolerance);
fprintf('RES 1 rows: L:%3i, W:%3i %s\n',L,W,fail);


 res = norm(fft(cast(f,'double'),2*L)-gga(f,linspace(0,1-0.5/L,2*L)));
 [test_failed,fail]=ltfatdiditfail(res,test_failed,tolerance);
 fprintf('RES 1/2 cols: L:%3i, W:%3i %s\n',L,W,fail);
 
 res = norm(fft(cast(f,'double'),5*L)-gga(f,linspace(0,1-0.2/L,5*L)));
 [test_failed,fail]=ltfatdiditfail(res,test_failed,tolerance);
 fprintf('RES 1/5 cols: L:%3i, W:%3i %s\n',L,W,fail);

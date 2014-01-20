function test_failed=test_signals;
%TEST_SIGNALS
%
%  The script ensures that all files are correctly included by loading
%  the singals and checking their sizes.
  
disp(' ===============  TEST_SIGNALS ===========');

test_failed=0;
  
[test_failed,fail]=ltfatdiditfail(numel(bat)-400,test_failed);
disp(['SIGNALS BAT       ',fail]);

[test_failed,fail]=ltfatdiditfail(numel(batmask)-1600,test_failed);
disp(['SIGNALS BATMASK   ',fail]);

[test_failed,fail]=ltfatdiditfail(numel(greasy)-5880,test_failed);
disp(['SIGNALS GREASY    ',fail]);

[test_failed,fail]=ltfatdiditfail(numel(linus)-41461,test_failed);
disp(['SIGNALS LINUS     ',fail]);

[test_failed,fail]=ltfatdiditfail(numel(gspi)-262144,test_failed);
disp(['SIGNALS GSPI      ',fail]);

[test_failed,fail]=ltfatdiditfail(numel(traindoppler)-157058,test_failed);
disp(['SIGNALS TRAINDOPPLER',fail]);

[test_failed,fail]=ltfatdiditfail(numel(otoclick)-2210,test_failed);
disp(['SIGNALS OTOCLICK  ',fail]);

[test_failed,fail]=ltfatdiditfail(numel(cameraman)-65536,test_failed);
disp(['SIGNALS CAMERAMAN ',fail]);

[test_failed,fail]=ltfatdiditfail(numel(cocktailparty)-363200,test_failed);
disp(['SIGNALS COCKTAILPARTY ',fail]);

[test_failed,fail]=ltfatdiditfail(numel(lichtenstein)-262144*3,test_failed);
disp(['SIGNALS LICHTENSTEIN ',fail]);


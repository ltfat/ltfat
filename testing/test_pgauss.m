function test_failed=test_pgauss
%TEST_PGAUSS  Test PGAUSS
  
    test_failed=0;
  
    disp(' ===============  TEST_PGAUSS ================');
    
    L=19;
    
    % Test that tfr=1 works
    res=norm(pgauss(L)-dft(pgauss(L)));
    [test_failed,fail]=ltfatdiditfail(res,test_failed);
    
    fprintf(['PGAUSS 1 %0.5g %s\n'],res,fail);
    
    % Test dilation property
    res=norm(pgauss(L,7)-dft(pgauss(L,1/7)));
    [test_failed,fail]=ltfatdiditfail(res,test_failed);
    
    fprintf(['PGAUSS 2 %0.5g %s\n'],res,fail);
    
    % Test norm
    res=norm(pgauss(L))-1;
    [test_failed,fail]=ltfatdiditfail(res,test_failed);
    
    fprintf(['PGAUSS 3 %0.5g %s\n'],res,fail);

    
    % Test that dft(freq shift) == time shift
    res=norm(dft(pgauss(L,'cf',5))-pgauss(L,'delay',5));
    [test_failed,fail]=ltfatdiditfail(res,test_failed);
    
    fprintf(['PGAUSS 3 %0.5g %s\n'],res,fail);

    


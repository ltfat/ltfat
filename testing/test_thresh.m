function test_failed=test_thresh
%TEST_THRESH  Compare sparse and full thesholding

test_failed=0;
disp(' ===============  TEST_THRESH ================');
global LTFAT_TEST_TYPE;
if ~strcmpi(LTFAT_TEST_TYPE,'double')
   disp(sprintf('Skipping. Cannot work with sparse matrices of type %s.',LTFAT_TEST_TYPE));
   return;
end

lambda=0.1;

ttypes={'hard','soft','wiener'};

for ii=1:2
  
  if ii==1
    g=tester_rand(3,4);
    field='REAL  ';
    g(2,2)=lambda;
  else
    g=tester_crand(3,4);
    field='CMPLX ';
    g(2,2)=lambda;
  end;
  
  for jj=1:3
    ttype=ttypes{jj};
    
    [xo_full, Nfull] = thresh(g,lambda,ttype,'full');
    [xo_sparse, Nsp] = thresh(g,lambda,ttype,'sparse');
    
    res = xo_full-xo_sparse;
    res = norm(res(:));
    
    res2 = Nfull-Nsp;
    
    [test_failed,fail]=ltfatdiditfail(res,test_failed);
    s=sprintf(['THRESH   %s %s %0.5g %s'],field,ttype,res,fail);
    disp(s);      
    
    [test_failed,fail]=ltfatdiditfail(res2,test_failed);
    s=sprintf(['THRESH N %s %s %0.5g %s'],field,ttype,res2,fail);      
    disp(s);
    
  end;
  
end;


% TEST_THRESH  

test_failed=0;

disp(' ===============  TEST_THRESH ================');

lambda=0.1;

for ii=1:2
  
  if ii==1
    g=rand(3,4);
    field='REAL  ';
    g(2,2)=lambda;
  else
    g=crand(3,4);
    field='CMPLX ';
    g(2,2)=lambda;
  end;
  
  for jj=1:2
    if jj==1
      ttype='hard';
    else
      ttype='soft';
    end;
    
    for kk=1:2

      if kk==1
        mtype='full';
        issparseval=0;
      else
        mtype='sparse';
        issparseval=1;
      end;
    
      xo     =     thresh(ttype,g,lambda,mtype);
      xo_ref = ref_thresh(ttype,g,lambda,mtype);
      
      res = xo-xo_ref;
      res = norm(res(:));
            
      [test_failed,fail]=ltfatdiditfail(res,test_failed);
      s=sprintf(['REF         %s %s %s %0.5g %s'],field,ttype,mtype,res,fail);
      disp(s);      

      res=(issparse(xo)==issparseval)-1;
      [test_failed,fail]=ltfatdiditfail(res,test_failed);
      s=sprintf(['REF SPARSE  %s %s %s %0.5g %s'],field,ttype,mtype,res,fail);
      disp(s);            
      
      [xo,N]         =     thresh(ttype,g,lambda,mtype);
      [xo_ref,N_ref] = ref_thresh(ttype,g,lambda,mtype);
      
      res1 = xo-xo_ref;
      res1 = norm(res1(:));

      res2 = N-N_ref;
      
      [test_failed,fail]=ltfatdiditfail(res1,test_failed);
      s=sprintf(['REF WITH N  %s %s %s %0.5g %s'],field,ttype,mtype,res1,fail);      
      disp(s);

      [test_failed,fail]=ltfatdiditfail(res2,test_failed);
      s=sprintf(['REF N VALUE %s %s %s %0.5g %s'],field,ttype,mtype,res2,fail);      
      disp(s);

      res=(issparse(xo)==issparseval)-1;      
      [test_failed,fail]=ltfatdiditfail(res,test_failed);
      s=sprintf(['REF SPARSE  %s %s %s %0.5g %s'],field,ttype,mtype,res,fail);
      disp(s);            

      
    end
    
  end;
  
end;


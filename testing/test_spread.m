function test_failed=test_spread
%TEST_SPREAD  Test spreading function.
%
%  This script runs a throrough test of the SPREADOP routine,
%  testing it on a range of input parameters.
global LTFAT_TEST_TYPE;
disp(' ===============  TEST_SPREAD ================');

disp('--- Used subroutines ---');

which comp_col2diag

Lr=[12, 13 ];

test_failed=0;

condNoLim.('single') = 100;
condNoLim.('double') = 100;

for ii=1:length(Lr);

L=Lr(ii); 

spfraction=.5;

for rtype=1:2
  
  if rtype==1
    rname='REAL ';	
  else
    rname='CMPLX';	
  end;
  
  for sptype=1:2
    
    if sptype==1
      spname='FULL  ';	
    else
      spname='SPARSE';	
      if ~strcmpi(LTFAT_TEST_TYPE,'double')
          disp(sprintf('Skipping. Cannot work with sparse matrices of type %s.',LTFAT_TEST_TYPE));
          break;
      end
    end;

    
    condNo = 1e10;
    while condNo>condNoLim.(LTFAT_TEST_TYPE)
    if rtype==1
      if sptype==1
        coef=tester_rand(L,L);
        T=tester_rand(L,L);
        coef2=tester_rand(L,L);
      else
        coef=tester_sprand(L,L,spfraction);
        T=tester_sprand(L,L,spfraction);
        coef2=tester_sprand(L,L,spfraction);          
      end;
    else
      if sptype==1
        coef=tester_crand(L,L);
        T=tester_crand(L,L);
        coef2=tester_crand(L,L);
      else
        coef=tester_sprand(L,L,spfraction);
        T=tester_sprand(L,L,spfraction);
        coef2=tester_sprand(L,L,spfraction);          
      end;
    end;
    
      coeftmp=ifft(full(coef))*L;
      % The following matrix is inverted in spreadinv. We want a nicer cond
      % number in order not to fail
      Ttmp=comp_col2diag(coeftmp);
      condNo = cond(Ttmp); 
    end

    
    for W=1:3
      
      if rtype==1
        f=tester_rand(L,W);
      else
        f=tester_crand(L,W);
      end;
      
      
      % ---------- Reference testing ------------ 
      
      fs=spreadop(f,coef);  
      fs2=ref_spreadop(f,coef,1);
      
      fsdiff=fs-fs2;
      res=norm(fsdiff(:));      
      [test_failed,fail]=ltfatdiditfail(res,test_failed);
      s=sprintf('REF %s %s L:%3i W:%2i %0.5g %s',rname,spname,L,W,res,fail);
      disp(s)

      % -------- Inversion testing -------------
     
      r=spreadinv(fs,coef);

      rdiff=f-r;
      res=norm(rdiff(:));      
      [test_failed,fail]=ltfatdiditfail(res,test_failed);
      s=sprintf('INV %s %s L:%3i W:%2i %0.5g %s',rname,spname,L,W,res,fail);
      disp(s)
      

      % -------- Twisted convolution -------------------

      f1=spreadop(spreadop(f,coef2),coef);
      f2=spreadop(f,tconv(coef,coef2));
      
      rdiff=f1-f2;
      res=norm(rdiff(:));      
      [test_failed,fail]=ltfatdiditfail(res,test_failed);
      s=sprintf('TWI %s %s L:%3i W:%2i %0.5g %s',rname,spname,L,W,res,fail);
      disp(s)
      
      % -------- Spreading function ---------------------
      
      coef=spreadfun(T);
      
      rdiff=T*f-spreadop(f,coef);
      
      res=norm(rdiff(:));      
      [test_failed,fail]=ltfatdiditfail(res,test_failed);
      s=sprintf('FUN %s %s L:%3i W:%2i %0.5g %s',rname,spname,L,W,res,fail);
      disp(s)

      % -------- Adjoint operator -----------------------
      
      cadj=spreadadj(coef);

      rdiff=T'*f-spreadop(f,cadj);
      
      res=norm(rdiff(:));      
      [test_failed,fail]=ltfatdiditfail(res,test_failed);
      s=sprintf('ADJ %s %s L:%3i W:%2i %0.5g %s',rname,spname,L,W,res,fail);
      disp(s)

      
    end;

  end;  
  
end;
end;



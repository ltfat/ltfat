function test_failed=test_dgt_fb_alg
%TEST_DGT_FB  Test the filter bank algorithms in DGT
%
%  This script runs a throrough test of the DGT routine,
%  testing it on a range of input parameters.
%
%  The script test the filter bank algorithms in DGT, IDGT, GABDUAL and
%  GABTIGHT by comparing with the full window case.
%
%  The computational backend is tested this way, but the
%  interfaces is not.
      
Lr  = [24, 35, 35, 24,144,108,144,135,77,77];
ar  = [ 6,  5,  5,  4,  9,  9, 12,  9, 7, 7];
Mr  = [ 8,  7,  7,  6, 16, 12, 24,  9,11,11];
glr = [16, 14, 21, 12, 48, 12, 24, 18,22,11];

test_failed=0;

disp(' ===============  TEST_DGT_FB_ALG ================');

disp('--- Used subroutines ---');

for ii=1:length(Lr);

  L=Lr(ii);
  
  M=Mr(ii);
  a=ar(ii);
  gl=glr(ii);

  b=L/M;
  N=L/a;
  

  for rtype=1:2
          
    if rtype==1
      rname='REAL ';	
      g=tester_rand(gl,1);
    else
      rname='CMPLX';	
      g=tester_crand(gl,1);
    end;
    
    if rtype==1
      rname='REAL ';	
      f=tester_rand(L,1);
    else
      rname='CMPLX';	
      f=tester_crand(L,1);
    end;
    
    cc = dgt(f,fir2long(g,L),a,M);
    
    cc_ref = ref_dgt_6(f,g,a,M);
    
    cdiff=cc-cc_ref;
    res=norm(cdiff(:));      

    [test_failed,fail]=ltfatdiditfail(res,test_failed);
    s=sprintf('REF  %s L:%3i a:%3i M:%3i gl:%3i %0.5g %s',rname,L,a,M,gl,res,fail);
    disp(s)
    
  end;

end;



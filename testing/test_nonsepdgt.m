function test_failed=test_nonsepdgt
%TEST_NONSEPDGT  Test non-separable DGT
%
%  This script runs a throrough test of the DGT routine,
%  testing it on a range of input parameters.
%
%  The computational backend is tested this way, but the
%  interface is not.
%
%  The script tests dgt, idgt, gabdual and gabtight.
%
%  Use TEST_WFAC and TEST_DGT_FAC for more specific testing
%  of the DGT backend.
      
Lr=[24,24];
ar=[ 4, 4];
Mr=[ 6, 6];
s1=[ 0, 1];
s2=[ 1, 2];

test_failed=0;

disp(' ===============  TEST_DGT ================');

disp('--- Used subroutines ---');

which comp_wfac
which comp_iwfac
which comp_dgt_long
which comp_idgt_fac
which comp_dgtreal_long
which comp_idgtreal_fac
which comp_gabdual_long
which comp_gabtight_long


for ii=1:length(Lr);

  L=Lr(ii);
  
  M=Mr(ii);
  a=ar(ii);
  s=[s1(ii), s2(ii)];
  
  b=L/M;
  N=L/a;
  c=gcd(a,M);
  d=gcd(b,N);
  p=a/c;
  q=M/c;
  

  for rtype=1:2
      
    if rtype==1
      rname='REAL ';	
      g=rand(L,1);
    else
      rname='CMPLX';	
      g=crand(L,1);
    end;
    
    gd=gabdual(g,a,M);
    %gt=gabtight(g,a,M);
    

    for W=1:3
          
      if rtype==1
        f=rand(L,W);
      else
        f=crand(L,W);
      end;      
      
      cc = nonsepdgt(f,g,a,M,s);
      
      cc_ref = ref_nonsepdgt(f,g,a,M,s);
      
      res = norm(cc(:)-cc_ref(:));
      [test_failed,fail]=ltfatdiditfail(res,test_failed);
      stext=sprintf(['REF %s L:%3i W:%2i a:%3i M:%3i s1:%3i s2:%3i %0.5g %s'],...
                    rname,L,W,a,M,s(1),s(2),res,fail);
      disp(stext)
      
      
      r=inonsepdgt(cc,gd,a,s);  
      res=norm(f-r,'fro');
      
      [test_failed,fail]=ltfatdiditfail(res,test_failed);
      stext=sprintf(['REC %s L:%3i W:%2i a:%3i M:%3i s1:%3i s2:%3i %0.5g %s'],...
                rname,L,W,a,M,s(1),s(2),res,fail);
      disp(stext)
      
      
    end;

  end;  

end;


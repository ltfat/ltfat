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

Lr =[24,24,36,36,48];
ar =[ 4, 4, 4, 4, 4];
Mr =[ 6, 6, 6, 6, 6];
lt1=[ 0, 1, 1, 2, 1];
lt2=[ 1, 2, 3, 3, 4];
    
test_failed=0;

disp(' ===============  TEST_NONSEPDGT ================');

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
  lt=[lt1(ii), lt2(ii)];
  
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
    
    gd=nonsepgabdual(g,a,M,lt);
    gd_smith=nonsepgabdual(g,a,M,lt,'smith');
    gd_shear=nonsepgabdual(g,a,M,lt,'shear');

    gt=nonsepgabtight(g,a,M,lt);
    gt_smith=nonsepgabtight(g,a,M,lt,'smith');
    gt_shear=nonsepgabtight(g,a,M,lt,'shear');

    for W=1:1
          
      if rtype==1
        f=rand(L,W);
      else
        f=crand(L,W);
      end;      
      
      % --------- test reference comparison ------------
      
      cc = nonsepdgt(f,g,a,M,lt);
      
      cc_ref = ref_nonsepdgt(f,g,a,M,lt);
      
      res = norm(cc(:)-cc_ref(:))/norm(cc(:));
      [test_failed,fail]=ltfatdiditfail(res,test_failed);
      stext=sprintf(['REF   %s L:%3i W:%2i a:%3i M:%3i lt1:%2i lt2:%2i %0.5g ' ...
      '%s'], rname,L,W,a,M,lt(1),lt(2),res,fail);
      disp(stext)

      
      % --------- test shear DGT -------------------------------
      
      cc_shear = nonsepdgt(f,g,a,M,lt,'shear');
            
      res = norm(cc(:)-cc_shear(:))/norm(cc(:));
      [test_failed,fail]=ltfatdiditfail(res,test_failed);
      stext=sprintf(['DGT SHREAR   %s L:%3i W:%2i a:%3i M:%3i lt1:%2i lt2:%2i %0.5g ' ...
      '%s'], rname,L,W,a,M,lt(1),lt(2),res,fail);
      disp(stext)

      % -------- test reconstruction using canonical dual -------
      
      r=inonsepdgt(cc,gd,a,lt);  
      res=norm(f-r,'fro');
      
      [test_failed,fail]=ltfatdiditfail(res,test_failed);
      stext=sprintf(['REC D %s L:%3i W:%2i a:%3i M:%3i lt1:%2i lt2:%2i %0.5g ' ...
      '%s'], rname,L,W,a,M,lt(1),lt(2),res,fail);
      disp(stext)
      
      % -------- test reconstruction using canonical dual -------

      cc = nonsepdgt(f,gt,a,M,lt);
      r=inonsepdgt(cc,gt,a,lt);  
      res=norm(f-r,'fro');
      
      [test_failed,fail]=ltfatdiditfail(res,test_failed);
      stext=sprintf(['REC T %s L:%3i W:%2i a:%3i M:%3i lt1:%2i lt2:%2i %0.5g ' ...
      '%s'], rname,L,W,a,M,lt(1),lt(2),res,fail);
      disp(stext)

      
    end;

    % -------- test frame bounds for tight frame -------
    
    B=nonsepgabframebounds(gt,a,M,lt);
    res=B-1;
    [test_failed,fail]=ltfatdiditfail(res,test_failed);
    stext=sprintf(['FRB   %s L:%3i a:%3i M:%3i lt1:%2i lt2:%2i %0.5g ' ...
                   '%s'], rname,L,a,M,lt(1),lt(2),res,fail);
    disp(stext)

    % -------- test smith and shear duals -------
    
    res=norm(gd-gd_smith)/norm(g);
    [test_failed,fail]=ltfatdiditfail(res,test_failed);
    stext=sprintf(['DUAL SMITH  %s L:%3i a:%3i M:%3i lt1:%2i lt2:%2i %0.5g ' ...
                   '%s'], rname,L,a,M,lt(1),lt(2),res,fail);
    disp(stext)
    res=norm(gd-gd_shear)/norm(g);
    [test_failed,fail]=ltfatdiditfail(res,test_failed);
    stext=sprintf(['DUAL SHEAR  %s L:%3i a:%3i M:%3i lt1:%2i lt2:%2i %0.5g ' ...
                   '%s'], rname,L,a,M,lt(1),lt(2),res,fail);
    disp(stext)
    
    % -------- test smith and shear tights -------
    
    res=norm(gt-gt_smith)/norm(g);
    [test_failed,fail]=ltfatdiditfail(res,test_failed);
    stext=sprintf(['TIGHT SMITH %s L:%3i a:%3i M:%3i lt1:%2i lt2:%2i %0.5g ' ...
                   '%s'], rname,L,a,M,lt(1),lt(2),res,fail);
    disp(stext)
    res=norm(gt-gt_shear)/norm(g);
    [test_failed,fail]=ltfatdiditfail(res,test_failed);
    stext=sprintf(['TIGHT SHEAR %s L:%3i a:%3i M:%3i lt1:%2i lt2:%2i %0.5g ' ...
                   '%s'], rname,L,a,M,lt(1),lt(2),res,fail);
    disp(stext)

  end;  

end;


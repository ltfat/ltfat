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

which comp_nonsepdgt_multi

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
  
  for gtype=1:2
      if gtype==1
          Lw=L;
      else
          Lw=M;
      end;
      
      for rtype=1:2
          
          if rtype==1
              rname='REAL';
              g=rand(Lw,1);
          else
              rname='CMPLX';	
              g=crand(Lw,1);
          end;
                    
          gd=gabdual(g,a,M,[],lt);
          gd_multi=gabdual(g,a,M,[],lt,'nsalg',1);
          gd_shear=gabdual(g,a,M,[],lt,'nsalg',2);

          gt=gabtight(g,a,M,[],lt);
          gt_multi=gabtight(g,a,M,[],lt,'nsalg',1);
          gt_shear=gabtight(g,a,M,[],lt,'nsalg',2);

          % For testing, we need to call some computational subroutines directly.
          gsafe=fir2long(g,L);
          gdsafe=fir2long(gd,L);

          
          for W=1:3
              
              if rtype==1
                  f=rand(L,W);
              else
                  f=crand(L,W);
              end;      
              
              % --------- test reference comparison ------------
              
              cc = dgt(f,g,a,M,[],lt);
              
              cc_ref = ref_nonsepdgt(f,g,a,M,lt);
              
              res = norm(cc(:)-cc_ref(:))/norm(cc(:));
              [test_failed,fail]=ltfatdiditfail(res,test_failed);
              stext=sprintf(['REF   %s L:%3i W:%2i LW:%3i a:%3i M:%3i lt1:%2i lt2:%2i %0.5g ' ...
                             '%s'], rname,L,W,Lw,a,M,lt(1),lt(2),res,fail);
              disp(stext)
              
              % --------- test multiwindow ---------------------
              
              cc = comp_nonsepdgt(f,gsafe,a,M,lt,0,1);
              
              res = norm(cc(:)-cc_ref(:))/norm(cc(:));
              [test_failed,fail]=ltfatdiditfail(res,test_failed);
              stext=sprintf(['DGT MULTIW %s L:%3i W:%2i LW:%3i a:%3i M:%3i lt1:%2i lt2:%2i %0.5g ' ...
                             '%s'], rname,L,W,Lw,a,M,lt(1),lt(2),res,fail);
              disp(stext)
              
              
              % --------- test shear DGT -------------------------------
              
              cc = comp_nonsepdgt(f,gsafe,a,M,lt,0,2);
              
              res = norm(cc(:)-cc_ref(:))/norm(cc(:));
              [test_failed,fail]=ltfatdiditfail(res,test_failed);
              stext=sprintf(['DGT SHREAR   %s L:%3i W:%2i LW:%3i a:%3i M:%3i lt1:%2i lt2:%2i %0.5g ' ...
                             '%s'], rname,L,W,Lw,a,M,lt(1),lt(2),res,fail);
              disp(stext)
              
              % -------- test reconstruction using canonical dual -------
              
              r=idgt(cc,gd,a,[],lt);
              res=norm(f-r,'fro');
              
              [test_failed,fail]=ltfatdiditfail(res,test_failed);
              stext=sprintf(['REC D %s L:%3i W:%2i LW:%3i a:%3i M:%3i lt1:%2i lt2:%2i %0.5g ' ...
                             '%s'], rname,L,W,Lw,a,M,lt(1),lt(2),res,fail);
              disp(stext)
              
              % -------- test reconstruction using canonical dual, multiwin algorithm -------
              
              r=comp_inonsepdgt(cc,gdsafe,a,lt,0,1);  
              res=norm(f-r,'fro');
              
              [test_failed,fail]=ltfatdiditfail(res,test_failed);
              stext=sprintf(['REC MULTIW D %s L:%3i W:%2i LW:%3i a:%3i M:%3i lt1:%2i lt2:%2i %0.5g ' ...
                             '%s'], rname,L,W,Lw,a,M,lt(1),lt(2),res,fail);
              disp(stext)
              
              % -------- test reconstruction using canonical dual, shear algorithm -------
              
              r=comp_inonsepdgt(cc,gdsafe,a,lt,0,2);  
              res=norm(f-r,'fro');
              
              [test_failed,fail]=ltfatdiditfail(res,test_failed);
              stext=sprintf(['REC SHEAR D %s L:%3i W:%2i LW:%3i a:%3i M:%3i lt1:%2i lt2:%2i %0.5g ' ...
                             '%s'], rname,L,W,Lw,a,M,lt(1),lt(2),res,fail);
              disp(stext)
              
              
              % -------- test reconstruction using canonical tight -------
              
              cc = dgt(f,gt,a,M,[],lt);
              r=idgt(cc,gt,a,[],lt);  
              res=norm(f-r,'fro');
              
              [test_failed,fail]=ltfatdiditfail(res,test_failed);
              stext=sprintf(['REC T %s L:%3i W:%2i LW:%3i a:%3i M:%3i lt1:%2i lt2:%2i %0.5g ' ...
                             '%s'], rname,L,W,Lw,a,M,lt(1),lt(2),res,fail);
              disp(stext)
              
              
          end;
          
          % -------- test frame bounds for tight frame -------
          
          B=gabframebounds(gt,a,M,[],lt);
          res=B-1;
          [test_failed,fail]=ltfatdiditfail(res,test_failed);
          stext=sprintf(['FRB   %s L:%3i LW:%3i a:%3i M:%3i lt1:%2i lt2:%2i %0.5g ' ...
                         '%s'], rname,L,Lw,a,M,lt(1),lt(2),res,fail);
          disp(stext)
          
          % -------- test multiwin dual -------
          
          res=norm(gd-gd_multi)/norm(g);
          [test_failed,fail]=ltfatdiditfail(res,test_failed);
          stext=sprintf(['DUAL MULTI  %s L:%3i LW:%3i a:%3i M:%3i lt1:%2i lt2:%2i %0.5g ' ...
                         '%s'], rname,L,Lw,a,M,lt(1),lt(2),res,fail);
          disp(stext)
          
          % -------- test shear dual -------
          
          res=norm(gd-gd_shear)/norm(g);
          [test_failed,fail]=ltfatdiditfail(res,test_failed);
          stext=sprintf(['DUAL SHEAR  %s L:%3i LW:%3i a:%3i M:%3i lt1:%2i lt2:%2i %0.5g ' ...
                         '%s'], rname,L,Lw,a,M,lt(1),lt(2),res,fail);
          disp(stext)
          
          % -------- test shear tight -------
          
          res=norm(gt-gt_multi)/norm(g);
          [test_failed,fail]=ltfatdiditfail(res,test_failed);
          stext=sprintf(['TIGHT MULTI %s L:%3i LW:%3i a:%3i M:%3i lt1:%2i lt2:%2i %0.5g ' ...
                         '%s'], rname,L,Lw,a,M,lt(1),lt(2),res,fail);
          disp(stext)

          % -------- test shear tight -------
          
          res=norm(gt-gt_shear)/norm(g);
          [test_failed,fail]=ltfatdiditfail(res,test_failed);
          stext=sprintf(['TIGHT SHEAR %s L:%3i LW:%3i a:%3i M:%3i lt1:%2i lt2:%2i %0.5g ' ...
                         '%s'], rname,L,Lw,a,M,lt(1),lt(2),res,fail);
          disp(stext)

          % ---- Test gabdualnorm --------------------------------------
          res=gabdualnorm(g,gd,a,M,[],lt);
          [test_failed,fail]=ltfatdiditfail(res,test_failed);
          stext=sprintf(['DUALNORM1 %s L:%3i LW:%3i a:%3i M:%3i lt1:%2i lt2:%2i %0.5g ' ...
                         '%s'], rname,L,Lw,a,M,lt(1),lt(2),res,fail);
          
          disp(stext);
          
          [o1,o2]=gabdualnorm(g,gd,a,M,[],lt);
          res=o1-1+o2;
          [test_failed,fail]=ltfatdiditfail(res,test_failed);
          stext=sprintf(['DUALNORM2 %s L:%3i LW:%3i a:%3i M:%3i lt1:%2i lt2:%2i %0.5g ' ...
                         '%s'], rname,L,Lw,a,M,lt(1),lt(2),res,fail);

          disp(stext);

          
      end;  
      
  end;

end;


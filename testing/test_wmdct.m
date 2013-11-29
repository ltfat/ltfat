function test_failed=test_wmdct
% Test the algorithm using LONG windows.

which comp_dwiltiii
which comp_idwiltiii

disp(' ===============  TEST_WMDCT ================');

Lr=[4, 6, 8,12,16,12,18,32,30];
Mr=[2, 3, 2, 3, 4, 2, 3, 4, 3];

test_failed=0;

for ii=1:length(Lr);
  for W=1:3
    for ftype=1:2
      for wtype=1:2
	L=Lr(ii);
	M=Mr(ii);
	
	a=M;
      
	if wtype==1
	  % Full length window
	  g=pgauss(L);
	  gd=wildual(g,M);
          wtype='LONG';
	else
	  g=firwin('sqrthann',2*M,'2');
	  gd=g;
          wtype='FIR ';
	end;
	
	if ftype==1
	  % Complex-valued test case
	  f=tester_crand(L,W);
	  S='CMPLX';
	else
	  % Real-valued tes
	  f=tester_rand(L,W);
	  S='REAL ';
	end;
	
	c=wmdct(f,g,M,L);  
	
	a=M;
	
	c2=ref_dwiltiii(f,g,a,M);
	r=iwmdct(c,gd);  
	
	res=norm(c(:)-c2(:));
	
        [test_failed,fail]=ltfatdiditfail(res,test_failed);        
	s=sprintf('REF  %s %s L:%3i W:%2i a:%3i M:%3i %0.5g %s',S,wtype,L,W,a,M,res,fail);
	disp(s);
	
	rdiff=f-r;
	res=norm(rdiff(:));
	
        [test_failed,fail]=ltfatdiditfail(res,test_failed);
        
	s=sprintf('REC  %s %s L:%3i W:%2i a:%3i M:%3i %0.5g %s',S,wtype,L,W,a,M,res,fail);
	disp(s);
        
        g=wilorth(M,L);
        c=wmdct(f,g,M);  
        r=iwmdct(c,g);
        rdiff=f-r;
        
	res=norm(rdiff(:));
	
        [test_failed,fail]=ltfatdiditfail(res,test_failed);
	s=sprintf('ORTH %s %s L:%3i W:%2i a:%3i M:%3i %0.5g %s',S,wtype,L,W,a,M,res,fail);
	disp(s);

	
      end;
    end;
  end;
end;











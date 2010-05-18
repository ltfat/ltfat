function test_failed=test_dwilt

Lr=[4, 6, 8,12,16,12,18,32,30];
Mr=[2, 3, 2, 3, 4, 2, 3, 4, 3];

disp(' ===============  TEST_DWILT ================');

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
	  g=firwin('sqrthann',2*M);
	  gd=g;
          wtype='FIR';
	end;
	
	if ftype==1
	  % Complex-valued test case
	  f=crand(L,W);
	  S='CMPLX';
	else
	  % Real-valued tes
	  f=rand(L,W);
	  S='REAL ';
	end;
	
	c=dwilt(f,g,M,L);  
      
	a=M;
	
	c2=ref_dwilt(f,g,a,M);
	r=idwilt(c,gd);  
	
	res=norm(c(:)-c2(:));
	
        fail='';
        if res>10e-10
	  fail='FAILED';
	  test_failed=test_failed+1;
        end;
        
	s=sprintf('REF  %s %s L:%3i W:%2i a:%3i M:%3i %0.5g %s',S,wtype,L,W,a,M,res,fail);
	disp(s)
	
	rdiff=f-r;
	res=norm(rdiff(:));
	
        fail='';
        if res>10e-10
	  fail='FAILED';
	  test_failed=test_failed+1;
        end;
   
    
	s=sprintf('REC  %s %s L:%3i W:%2i a:%3i M:%3i %0.5g %s',S,wtype,L,W,a,M,res,fail);
	disp(s)

	g=wilorth(M,L);
        c=dwilt(f,g,M);  
        r=idwilt(c,g);
        rdiff=f-r;
        
	res=norm(rdiff(:));
	
        fail='';
        if res>10e-10
	  fail='FAILED';
	  test_failed=test_failed+1;
        end;
       
	s=sprintf('ORTH %s %s L:%3i W:%2i a:%3i M:%3i %0.5g %s',S,wtype,L,W,a,M,res,fail);
	disp(s)
        
                        
      end;
    end;
  end;
end;











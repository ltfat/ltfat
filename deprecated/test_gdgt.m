function test_failed=test_gdgt
%TEST_GDGT  Test GDGT
%
%  This script runs a throrough test of the COMP_GDGT routine, testing it on
%  a range of input parameters.
%


Lr=[24,144,108,144,24,135,35,77,20];
ar=[ 4,  9,  9, 12, 6,  9, 5, 7, 1];
Mr=[ 6, 16, 12, 24, 8,  9, 7,11,20];

R=1;

test_failed=0;

disp(' ===============  TEST_GDGT ================');

for ii=1:length(Lr);

  for ctii=0:1

    c_t=ctii*.5;

    for cfii=0:1      
	
      c_f=cfii*.5;

      for W=1:3
	
	L=Lr(ii);
	
	M=Mr(ii);
	a=ar(ii);
	
	b=L/M;
	N=L/a;
	c=gcd(a,M);
	d=gcd(b,N);
	p=a/c;
	q=M/c;
	
	%g=(1:L)';
	%g=i*gabtight(a,b,L);
	f=tester_crand(L,W);
	g=tester_crand(L,R);
	
	gd=gabdual(g,a,M);
	gt=gabtight(g,a,M);
	
	cc=comp_gdgt(f,g,a,M,L,c_t,c_f,0,0);  
	cc2=reshape(ref_gdgt(f,g,a,M,c_t,c_f,0),M,N,W);
	
	cdiff=cc-cc2;        
	res=norm(cdiff(:));      
        [test_failed,fail]=ltfatdiditfail(res,test_failed);
	s=sprintf(['REF L:%3i c_t: %0.5g c_f %0.5g ', ...
                   'W:%2i R:%2i a:%3i b:%3i c:%3i ', ...
                   'd:%3i p:%3i q:%3i %0.5g %s'],...
                  L,c_t,c_f,W,R,a,b,c,d,p,q,res,fail);
	disp(s)
	
	r=comp_igdgt(cc,gd,a,M,L,c_t,c_f,0,0);  
	res=norm(f-r,'fro');
	[test_failed,fail]=ltfatdiditfail(res,test_failed);
        
	s=sprintf(['REC L:%3i c_t: %0.5g c_f %0.5g ',...
                   'W:%2i R:%2i a:%3i b:%3i c:%3i ',...
                   'd:%3i p:%3i q:%3i %0.5g %s'],...
                  L,c_t,c_f,W,R,a,b,c,d,p,q,res,fail);
	disp(s)
	
	res=norm(f-idgt(dgt(f,gt,a,M),gt,a),'fro');
        [test_failed,fail]=ltfatdiditfail(res,test_failed);
	s=sprintf(['TIG L:%3i c_t: %0.5g c_f %0.5g ',...
                   'W:%2i R:%2i a:%3i b:%3i c:%3i ',...
                   'd:%3i p:%3i q:%3i %0.5g %s'],...
                  L,c_t,c_f,W,R,a,b,c,d,p,q,res,fail);
	disp(s);
      end;
      
    end;
  end;
end;






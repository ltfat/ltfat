function test_failed=test_pfilt_1
%
%
%   This is the old test_pfilt from before the struct filters was
%   introduced.

Lr =[9,9,10,10,10,12];
Lgr=[9,4,10, 7,10,12];
ar =[3,3, 5, 5, 1, 3];

test_failed=0;

disp(' ===============  TEST_PFILT_1 ==============');

disp('--- Used subroutines ---');

which comp_pfilt

for jj=1:length(Lr)
  L=Lr(jj);
  Lg=Lgr(jj);
  a=ar(jj);
  
  for W=1:3
  
    for rtype=1:2
      if rtype==1
        rname='REAL ';	
        f=tester_rand(L,W);
        g=tester_rand(Lg,1);
      else
        rname='CMPLX';	
        f=tester_crand(L,W);
        g=tester_crand(Lg,1);
      end;
                 
      h1=pfilt(f,g,a);
      h2=ref_pfilt(f,g,a);
      
      res=norm(h1-h2);
      [test_failed,fail]=ltfatdiditfail(res,test_failed);        
      s=sprintf('PFILT %3s  L:%3i W:%3i Lg:%3i a:%3i %0.5g %s',rname,L,W,Lg,a,res,fail);
      disp(s);
    end;
  end;
end;



function test_failed=test_frametf
%TEST_FRAMETF  Test the frames tf-plane conversion
%
%   This tests if framecoef2tf and frametf2coef work.

test_failed=0;
L = 456;
W = 3;

f = tester_rand(L,W);

Fr{1}  = frame('dgt','gauss',10,20);
Fr{1}  = frame('dgtreal','gauss',10,20);
Fr{3}  = frame('dwilt','gauss',20);
Fr{4}  = frame('wmdct','gauss',20);

   gfilt={tester_rand(30,1),...
          tester_rand(20,1),...
          tester_rand(15,1),...
          tester_rand(10,1)};
      
Fr{5} = frame('ufilterbank',    gfilt,3,4);

Fr{6} = frame('ufwt','db4',4);

Fr{7} = frame('uwfbt',{'db4',4});

Fr{8} = frame('uwpfbt',{'db4',4});

for ii=1:numel(Fr)
  
  F=Fr{ii};
  
  % To avoid holes in Fr
  if isempty(F)
    continue;
  end;
  
  c = frana(F,f);
  ctf = framecoef2tf(F,c);
  
  c2 = frametf2coef(F,ctf);
  
  res = norm(c-c2);
  [test_failed,fail]=ltfatdiditfail(res,test_failed);
  fprintf('COEFEQ  %s  %0.5g %s\n',F.type,res,fail);
 
end

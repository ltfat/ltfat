function test_failed=test_gabmuleigs
%TEST_GABMULEIGS  Test GABMULEIGS
%
%   Test GABMULEIGS by comparing the output from the iterative and full algorithm.

disp(' ===============  TEST_GABMULEIGS ================');

test_failed=0;
  
a=20;
M=30;

L=a*M;
N=L/a;

c=randn(M,N);

g=gabtight(a,M,L);

% [V1,D1]=gabmuleigs(10,c,g,a,'iter');
% [V2,D2]=gabmuleigs(10,c,g,a,'full');
F = frame('dgt',g,a,M);
c = framenative2coef(F,c);
[V1,D1]=framemuleigs(F,F,c,10,'iter');
[V2,D2]=framemuleigs(F,F,c,10,'full');

res=norm(D1-D2);

[test_failed,fail]=ltfatdiditfail(res,test_failed);
s=sprintf('GABMULEIGS   L:%3i a:%3i M:%3i %0.5g %s',L,a,M,res,fail);
disp(s);





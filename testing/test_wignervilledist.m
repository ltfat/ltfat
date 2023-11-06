function test_failed = test_wignervilledist
Lr = [1, 19, 20];

test_failed = 0;

disp(' ===============  TEST_WIGNERVILLEDIST ==============');

for ii = 1: length(Lr)
  L = Lr(ii);
    for n = 1:4
    
    if (n==1)
    type1 = 'auto';
    type2 = 'real';
    f = tester_rand(L,1);
    g = f;
    elseif (n==2)
    type1 = 'auto';
    type2 = 'complex';
    f = tester_crand(L,1);
    g = f;
    elseif (n==3)
    type1 = 'cross';
    type2 = 'real';
    f = tester_rand(L,1);
    g = tester_rand(L,1);
    elseif (n==4)
    type1 = 'cross';
    type2 = 'complex';
    f = tester_crand(L,1);
    g = tester_crand(L,1);
    end
  
    r1 = ref_wignervilledist(f,g);
    r2 = wignervilledist(f,g);
  
    res = norm(r1-r2);
  
    [test_failed, fail] = ltfatdiditfail(res, test_failed);
    s = sprintf('DWVD %3s %3s L:%3i %0.5g %s', type1, type2, L, res, fail);
    disp(s);
    end
end
end

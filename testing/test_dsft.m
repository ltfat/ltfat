function test_failed = test_dsft
Lr = [2, 19, 20];

test_failed = 0;

disp(' ===============  TEST_DSFT ==============');

for ii = 1: length(Lr)
  L = Lr(ii);
    for n = 1:4
    
    if (n == 1)
    type = 'real';
    L2 = L;
    F = tester_rand(L,L);
    elseif (n == 2)
    type = 'real';
    L2 = L+1;
    F = tester_rand(L,L2);
    elseif (n == 3)
    type = 'complex';
    L2 = L;
    F = tester_crand(L,L);
    else
    type = 'complex';
    L2 = L+1;
    F = tester_crand(L,L2); 
    end
    
    r1 = ref_dsft(F);
    r2 = dsft(F);
  
    res = norm(r1-r2);
  
    [test_failed, fail] = ltfatdiditfail(res, test_failed);
    s = sprintf('DSFT %3s size: %dx%d %0.5g %s', type, L, L2, res, fail);
    disp(s);
    end
end

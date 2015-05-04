function test_failed = test_drihaczekdist
Lr = [1, 19, 20];

test_failed = 0;

disp(' ===============  TEST_DRIHACZEKDIST ==============');

for ii = 1: length(Lr)
  L = Lr(ii);
    for type = {'real', 'complex'}
    
    if strcmp(type{1}, 'real')
    f = tester_rand(L,1);
    else
    f = tester_crand(L,1);
    end
  
    r1 = ref_drihaczekdist(f);
    r2 = drihaczekdist(f);
  
    res = norm(r1-r2);
  
    [test_failed, fail] = ltfatdiditfail(res, test_failed);
    s = sprintf('DRIHACZEKDIST %3s L:%3i %0.5g %s', type{1}, L, res, fail);
    disp(s);
    end
end

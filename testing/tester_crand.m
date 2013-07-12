function f=tester_crand(p1,p2);
%CRAND   Random complex numbers for testing.
%   Usage: f=tester_crand(p1,p2);

global LTFAT_TEST_TYPE

if isempty(LTFAT_TEST_TYPE)
    LTFAT_TEST_TYPE='double';
end;

f=rand(p1,p2)-.5+i*(rand(p1,p2)-.5);

if strcmp(LTFAT_TEST_TYPE,'single')
    f=single(f);
end;
    


function f=tester_rand(varargin);
%RAND   Random numbers for testing.
%   Usage: f=tester_rand(p1,p2);

global LTFAT_TEST_TYPE

if isempty(LTFAT_TEST_TYPE)
    LTFAT_TEST_TYPE='double';
end;

f=rand(varargin{:}) - 0.5;

if strcmp(LTFAT_TEST_TYPE,'single')
    f=single(f);
end;
    


function f=tester_sprand(varargin);
%TESTER_SPRAND   Sparse random numbers for testing.
%   Usage: f=tester_sprand(p1,p2);

global LTFAT_TEST_TYPE

if isempty(LTFAT_TEST_TYPE)
    LTFAT_TEST_TYPE='double'
end;

f=sprand(varargin{:});

if strcmp(LTFAT_TEST_TYPE,'single')
    f=single(f);
end;
    


function test_failed = test_gabwintight()
%TEST_GABWINTIGHT Some of the windows returned from firwin are tight
%   immediatelly. This function tests for that
test_failed = 0;



disp(' ===============  TEST_GABWINTIGHT ================');

disp('--- Used subroutines ---');

which gabwin
which comp_window


shouldBeTight = [...
    struct('g','sine','a',10,'M',40,'L',[]),...
    struct('g',{{'sine',40}},'a',10,'M',40,'L',[]),...
    struct('g',{{'sine',28}},'a',14,'M',40,'L',[]),...
    struct('g','sine','a',20,'M',40,'L',[]),...
    struct('g',{{'sine',60}},'a',20,'M',60,'L',[]),...
    struct('g',{{'sine',54}},'a',18,'M',60,'L',[]),...
    struct('g',{{'sine',54,'inf'}},'a',18,'M',60,'L',[]),...
    struct('g','sine','a',5,'M',40,'L',[]),...
    struct('g','sqrttria','a',5,'M',40,'L',[]),...
];

shouldNotBeTight = [...
    struct('g','sine','a',16,'M',40,'L',[]),...
    struct('g',{{'sine',41}},'a',10,'M',40,'L',[]),...
    struct('g',{{'sine',26}},'a',10,'M',40,'L',[]),...
    struct('g','sqrttria','a',40,'M',40,'L',[]),...
];


for ii=1:numel(shouldBeTight)
    gw = shouldBeTight(ii);
    
    [~,info] = gabwin(gw.g,gw.a,gw.M,gw.L);
    [test_failed,fail]=ltfatdiditfail(~info.istight,test_failed,0);
    fprintf(['GABWINISTIGHT g= a=%i M=%i %s\n'],gw.a,gw.M,fail);
    
end

for ii=1:numel(shouldNotBeTight)
    gw = shouldNotBeTight(ii);
    
    [~,info] = gabwin(gw.g,gw.a,gw.M,gw.L);
    [test_failed,fail]=ltfatdiditfail(info.istight,test_failed,0);
    fprintf(['GABWINISNOTTIGHT g= a=%i M=%i %s\n'],gw.a,gw.M,fail);
    
end
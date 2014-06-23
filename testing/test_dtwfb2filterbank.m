function test_failed = test_dtwfb2filterbank
% This test only 
test_failed = 0;

L = 1000;


dualwt{1} = {'ana:qshift1',3};
dualwt{2} = {'ana:qshift3',5,'first','symorth1'};
dualwt{3} = {'ana:qshift3',5,'full','first','symorth1'};
dualwt{4} = {'ana:oddeven1',5};


for ii = 1:numel(dualwt)
    [g,a,info] = dtwfb2filterbank( dualwt{ii},'real');

    G = filterbankfreqz(g,a,L);

    Greal = filterbankfreqz(info.g1,a,L);
    Ghilb = filterbankfreqz(info.g2,a,L);

    G2 = Greal+1i*Ghilb;


    res = norm(G(:)-G2(:));

    [test_failed,fail]=ltfatdiditfail(res,test_failed);
    s=sprintf(['DUAL-TREE %i %0.5g %s'],ii,res,fail);
    disp(s)
end

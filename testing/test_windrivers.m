function test_failed=test_windrivers
%TEST_WINDRIVERS  Test if the window drivers pass certain construction
  
    test_failed=0;
  
    disp(' ===============  TEST_WINDRIVERS ================');
    
    
a=5;
M=6;
L=60;

% We expect that if the following commands finish, they produce the correct
% output, so we only test that they do not generate fatal errors.

g=gabwin('gauss',a,M,L);
g=gabwin({'gauss',1},a,M,L);
gd=gabwin('gaussdual',a,M,L);
gd=gabwin({'tight','gauss'},a,M,L);
g=gabwin({'dual',{'hann',M}},a,M,L);



function f=crand(p1,p2);
%CRAND   Random complex numbers for testing.
%   Usage: f=tester_crand(p1,p2);

f=rand(p1,p2)-.5+i*(rand(p1,p2)-.5);



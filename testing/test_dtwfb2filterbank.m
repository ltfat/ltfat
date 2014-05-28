function test_failed = test_dtwfb2filterbank
test_failed = 0;

L = 1000;

dualwt = {'qshift3',4,'first','db10'};

[g,a] = dtwfbreal2filterbank( dualwt,'nat');

G = filterbankfreqz(g,a,L);
G2 = calc_responses(dualwt, L,'nat');

res = norm(G(:)-G2(:));

[test_failed,fail]=ltfatdiditfail(res,test_failed);
s=sprintf(['DUAL-TREE  %0.5g %s'],res,fail);
disp(s)

figure(1);
plot(abs(G));

figure(2);
plot(abs(G2));

function G = calc_responses(dualwt, L,varargin)


dtw = dtwfbinit({'strict',dualwt},varargin{:});

wtPath = 1:numel(dtw.nodes);
wtPath(noOfNodeOutputs(1:numel(dtw.nodes),dtw)==0)=[];
rangeLoc = rangeInLocalOutputs(wtPath,dtw);
rangeOut = rangeInOutputs(wtPath,dtw);
[g,a] = nodesMultid(wtPath,rangeLoc,rangeOut,dtw);

dtw.nodes = dtw.dualnodes;
[g2,a2] = nodesMultid(wtPath,rangeLoc,rangeOut,dtw);


gf = filterbankfreqz(g,a,L);
gf2 = filterbankfreqz(g2,a2,L);

G = gf + i*gf2;



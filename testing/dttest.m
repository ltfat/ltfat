%f = randn(10000,1);
f = zeros(1500,1);
f(45) = 1;

wt =  {'qshift3',6,'first','ana:symorth3'};

[g,a,info] = dtwfbreal2filterbank( wt);
c = filterbank(f,g,a);

wt2 =  dtwfbinit(wt);

figure(1);
clf;
F =filterbankfreqz(g,a,2048,'linabs','plot');

%figure(2);clf;
%plot(pgrpdelay(g{1},2048));


cd = dtwfbreal(f,wt);
  

norm(cell2mat(c)-cell2mat(cd))
figure(4);
stem(abs([(cell2mat(c)),cell2mat(cd)]))

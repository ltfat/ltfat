%f = randn(2000,1);
f = zeros(2000,1);
f(149) = 1;


wt =  {'oddeven1',3};

wt2 =  dtwfbinit(wt);

%[g,a,info] = dtwfb2filterbank( wt,'real','nat');
%c = filterbank(f,g,a);



%figure(1);plot(pgrpdelay(wt2.nodes{2}.g{1},1024))
%figure(2);plot(pgrpdelay(wt2.dualnodes{2}.g{1},1024))
%figure(1);plot([pgrpdelay(info.g1{1},1024),pgrpdelay(info.g2{1},1024)]);
%figure(2);plot([pgrpdelay(info.g1{2},1024),pgrpdelay(info.g2{2},1024)]);


figure(1);
clf;
F =filterbankfreqz(g,a,2*2048,'linabs','plot');

[cd,info] = dtwfbreal(f,wt,'nat');
%figure(1);
%plotwavelets(cd,info,'dynrange',80);

%czero = cellfun(@(cEl) zeros(size(cEl)),cd,'UniformOutput',0);


fhat = idtwfbreal(cd,wt,numel(f),'nat');


fprintf('Rec. err.  %d\n',norm(f-fhat));
figure(5);
clf;
stem([f,fhat]);
  
fprintf('Coeff. err. %d\n',norm(cell2mat(c)-cell2mat(cd)));

figure(4);
stem(abs([(cell2mat(c)),cell2mat(cd)]))

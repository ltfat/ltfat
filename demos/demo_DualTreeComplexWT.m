%DEMO_DUALTREECOMPLEXWT Building Custom wavelet tree 

f = gspi;
J = 10;

% "Real" tree
realTree = wfbtinit({'dtree',1});
imagTree = wfbtinit({'dtree',2});

wother{1} =  fwtinit({'dtree',3});
wother{2} =  fwtinit({'dtree',4});

for jj=1:J-1
    realTree = wfbtput(jj,0,wother{2-rem(jj,2)},realTree);
    imagTree = wfbtput(jj,0,wother{1+rem(jj,2)},imagTree);
end


% figure(1);clf;
% freqzfb(multid(realTree));
% figure(2);clf;
% freqzfb(multid(imagTree));


creal = wfbt(f,realTree);
cimag = wfbt(f,imagTree);

ccmplx = cell(size(creal));

for jj=1:length(creal)
    ccmplx{jj} = creal{jj} + i*cimag{jj}; 
end

plotfwt(ccmplx);
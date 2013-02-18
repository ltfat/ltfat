function test_failed = test_wtfft_pr
%TEST_COMP_FWT_ALL
%
% Checks perfect reconstruction of the wavelet transform
%
test_failed = 0;



load vonkoch;
f=vonkoch;f=f';
%f = 0:6;
%f = flipud(f');
%f = ones(7,1);
%f = randn(32,1);

J = 7;

w = fwtinit({'dden',3});
%[w.h,w.g,abase]=wfilt_hden(1); 

c2 = fwt(f,w,J);
fhat2 = ifwt(c2,w,J,length(f));

[h,a] = wfbtmultid({w,J},'ana');
[g,a] = wfbtmultid({w,J},'syn');
figure(3);freqzfb(h,length(f));
figure(4);freqzfb(g,length(f));
H = freqzfb(h,filterbanklength(length(f),a));
G = freqzfb(g,filterbanklength(length(f),a));



c1 = wtfft(f,H,a);

figure(2);clf;plotwavc(c1);

fhat = iwtfft(c1,G,a,length(f));


figure(1);clf;stem([f,fhat,fhat2]);
legend({'orig','iwtfft'});
title(sprintf('norm(f-fhat)=%d',norm(f-fhat2)));

c2form = cell(numel(c2)-(length(w.h)-2),1);
c2form{1} = c2{1};
cSformIdx = 2;
for jj=2:J+1
    for ii=1:length(w.h)-1
       c2form{cSformIdx} = c2{jj,ii};
       cSformIdx=cSformIdx+1;
    end
end
figure(5);
printCoeffs( c1,c2form);


 
function printCoeffs( x,y)

[J,N1] = size(x);

for j=1:J
    subplot(J,1,j);
    % err = x{j}(:) - y{j}(:);
      stem([x{j}(:),y{j}(:)]);
      lh = line([0 length(x{j})],[eps eps]);
      set(lh,'Color',[1 0 0]);
      lh =line([0 length(x{j})],[-eps -eps]);
      set(lh,'Color',[1 0 0]);

end

function coefs = coefMatToLTFAT(C,S,lo_r,hi_r,J)

coefs = cell(J+1,1);

coefs{1,1} = appcoef(C,S,lo_r,hi_r,J);
for j=1:J
     [coefs{end-j+1}] = detcoef(C,S,j); 
end


 
   

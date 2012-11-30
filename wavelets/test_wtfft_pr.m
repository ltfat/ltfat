function test_failed = test_wtfft_pr
%TEST_COMP_FWT_ALL
%
% Checks perfect reconstruction of the wavelet transform
%
test_failed = 0;



load vonkoch;
f=vonkoch;
%f = 0:2^9-1;
f = f';

J = 2;


[w.h,w.g,abase]=wfilt_mband(1); 



[h,a] = multid(w.h,J,abase,'full');
[g,a] = multid(w.g,J,abase,'full','syn');
figure(3);freqzfb(h,length(f));
figure(4);freqzfb(g,length(f));
H = freqzfb(h,filterbanklength(length(f),a));
G = freqzfb(g,filterbanklength(length(f),a));


%[H,G] = wfreq_lemaire(length(f));

c1 = wtfft(f,H,a);

figure(2);clf;plotwavc(c1,'undec');

fhat = iwtfft(c1,G,a,length(f));


figure(1);clf;stem([f,fhat]);
legend({'orig','iwtfft'});
title(sprintf('norm(f-fhat)=%d',norm(f-fhat)));




 
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


 
   

function test_failed = test_wtfft
%TEST_COMP_FWT_ALL
%
% Checks perfect reconstruction of the wavelet transform
%
test_failed = 0;



load vonkoch;
f=vonkoch;
f = 0:2^9-1;
f = f';

J = 4;
wname = 'bior2.2';
%w = waveletfb('hden',1);
[hw,gw] = wfilt_matlabwtwrapper(wname);

dwtmode('per');
[C I] = wavedec(f,J,wname);
[lo,hi,lo_r,hi_r] = wfilters(wname);
cmat = coefMatToLTFAT(C,I,lo_r,hi_r,J)


[h,a] = multid(hw,J,[2,2]);
[g,a] = multid(gw,J,[2,2],'syn');
H = freqzfb(h,length(f));
G = freqzfb(g,length(f));
c1 = wtfft(f,H,a);
c2 = fwt(f,hw,J);

fhat = iwtfft(c1,G,a,length(f));
fhat2 = ifwt(c2,gw,length(f));

figure(1);clf;stem([f,fhat,fhat2]);legend({'orig','iwtfft','ifwt'});

c2form = cell(numel(c2),1);
c2form{1} = c2{1};
cSformIdx = 2;
for jj=2:J+1
    for ii=1:length(hw)-1
       c2form{cSformIdx} = c2{jj,ii};
       cSformIdx=cSformIdx+1;
    end
end

figure(2);clf;printCoeffs(c1,cmat);


 
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


 
   

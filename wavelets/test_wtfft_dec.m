function test_failed = test_wtfft_dec;
%TEST_COMP_FWT_ALL
%
% Checks perfect reconstruction of the wavelet transform
%
test_failed = 0;



load vonkoch;
f=vonkoch;
f = 0:2^8-1;
f = f';

J = 1;
w = waveletfb('algmband',2);


[h,a] = multid(w.h,J);
H = freqzfb(h,filterbanklength(length(f),a));
tic; c1 = wtfft(f,H,a);toc;
tic; c2 = fwt(f,w,J); toc;


c2form = cell(numel(c2)-(length(w.h)-2),1);
c2form{1} = c2{1};
cSformIdx = 2;
for jj=2:J+1
    for ii=1:length(w.h)-1
       c2form{cSformIdx} = c2{jj,ii};
       cSformIdx=cSformIdx+1;
    end
end
printCoeffs(c1,c2form);

% hlens = zeros(numel(h),1);
% for jj = 1:numel(h)
%     hlens(jj) = length(h{jj});
% end
% 
% shifts=zeros(numel(c1),1);
% % check coefficients and find
% for jj = 1:numel(c1)
%     if(norm(c1{jj}-c2{jj})>1e-6)
%         
%         for sh=1:floor(length(c1{jj})/2)
%            if(norm(c1{jj}-circshift(c2{jj},sh))<1e-6)
%              shifts(jj)=sh;
%              continue;
%            end 
%            
%            if(norm(c1{jj}-circshift(c2{jj},-sh))<1e-6)
%              shifts(jj)=-sh;
%              continue;
%            end 
%         end
%         
%         if(shifts(jj)~=0)  continue; end;
%         % even all coefficients shifts are not equal
%      shifts(jj)=Inf;
%     end
% end
% 
% shifts
% hlens
 
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






 
   

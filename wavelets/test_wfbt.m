function test_wfbt

% 'mband' 4 filters, subs 4
% w = fwtinit({'mband',1});

% 'mband' 4 filters, subs 2
% w = waveletfb({'dgrid',1});
 
% 'mband' 3 filters, subs 2,2,1
%w = waveletfb({'hden',3});

w = fwtinit({'db',10});


J = 5;
wt = wfbtinit(w,'J',J,'full');
wt = wfbtremove(2,0,wt,'force');
wt = nat2freqOrder(0,wt);
wtdual = nat2freqOrder(1,wt);


% Building custom filterbank tree
% DWT case
% Level J DWT
% for jj=1:J
%     wt = wfbtput(jj-1,0,w.h,w.a,wt);
% end
L = 1024;
f = randn(L,1);

f=gspi;
L = length(f);
figure(1);
freqzfb(multid(wtdual),1024*64);


%[h,a] = multid(wt);
%H = freqzfb(h,filterbanklength(length(f),a));

%c1 = wtfft(f,H,a);
tic;
c2 = wfbt(f,wt);
toc;

tic;
fhat = iwfbt(c2,wtdual,length(f));
toc;

errf = norm(f-fhat);
if(errf>1e-6)
  figure(2);
  clf;
  stem([f,fhat]);
end

% if(checkCoefs(c1,c2)>1e-10)
%     printCoeffs(c1,c2);
%   disp('Not OK');
% else
%     disp('OK');
% end
figure(1);
plotfwt(c2);


% %plot frequency response
% freqzfb(multid(wtree));
% 
% wtree = wfbtinit();
% % Full depth-J tree
% for jj=1:J
%     for kk=1:length(w.h)^(jj-1)
%        wtree = wfbtput(jj-1,kk-1,w.h,w.a,wtree);
%     end
% end
% 
% %plot frequency response
% freqzfb(multid(wtree));
% 
% % delete some nodes
% wtree = wfbtremove(1,1,wtree,'force');
% wtree = wfbtremove(J-1,0,wtree);
% 
% freqzfb(multid(wtree));
% 
% [h,a] = multid(wtree);
% H = freqzfb(h);

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


function err = checkCoefs(c,coefs)
err = 0;
for ii=1:length(c)
   err = err +  norm(c{ii}(:) - coefs{ii}(:));  
end








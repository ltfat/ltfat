function test_failed=test_undeceq
% This function test whether fwt is just a subsampled version of ufwt, wfbt
% of uwfbt etc.
test_failed = 0;

J = 5;

L = 128;

f = tester_rand(L,1);

wav = {'db4','spline4:4'};

for ii = 1:numel(wav)
    
   w = fwtinit(wav{ii});
    
   [c,info] = fwt(f,wav{ii},J,'cell');
   cu = ufwt(f,wav{ii},J,'noscale');
   
   err = 0;
   suFac = size(cu,1)./info.Lc;
   
   for jj = 1:numel(c)
       err = err + norm(c{jj}-cu(1:suFac(jj):end,jj));
   end
   
   [test_failed,fail]=ltfatdiditfail(err,test_failed);
    fprintf('DWT J=%d, %6.6s, L=%d, err=%.4e %s \n',J,wav{ii},length(f),err,fail);
end
   
for ii = 1:numel(wav)
   [c,info] = wfbt(f,{wav{ii},J});
   cu = uwfbt(f,{wav{ii},J},'noscale');
   
   err = 0;
   suFac = size(cu,1)./info.Lc;
   
   for jj = 1:numel(c)
       err = err + norm(c{jj}-cu(1:suFac(jj):end,jj));
   end
   
   [test_failed,fail]=ltfatdiditfail(err,test_failed);
    fprintf('WFBT J=%d, %6.6s, L=%d, err=%.4e %s \n',J,wav{ii},length(f),err,fail);
end

for ii = 1:numel(wav)
   [c,info] = wpfbt(f,{wav{ii},J});
   cu = uwpfbt(f,{wav{ii},J},'noscale');
   
   err = 0;
   suFac = size(cu,1)./info.Lc;
   
   for jj = 1:numel(c)
       err = err + norm(c{jj}-cu(1:suFac(jj):end,jj));
   end
   
   [test_failed,fail]=ltfatdiditfail(err,test_failed);
    fprintf('WPFBT J=%d, %6.6s, L=%d, err=%.4e %s \n',J,wav{ii},length(f),err,fail);
end

function test_failed = test_fwt2

disp('========= TEST FWT2 ============');
test_failed = 0;

dims = { [20,30], [150,151],[226,253], };
flags = {'standard','tensor'};
filt = {{'mband1',2},{'db10',4},{'sym8',4},{'spline4:4',4},};


for ii=1:numel(dims)
   f = randn(dims{ii});
   for jj=1:numel(flags)
      for ff=1:numel(filt)
         c = fwt2(f,filt{ff}{1},filt{ff}{2},flags{jj});
         fhat = ifwt2(c,filt{ff}{1},filt{ff}{2},dims{ii},flags{jj});
         err = norm(f-fhat,'fro');
         fprintf('J=%d, filt=%s, dim=[%d,%d], flag=%s, err=%d \n',filt{ff}{2},filt{ff}{1},size(f,1),size(f,2),flags{jj},err);
         if err>1e-6
            test_failed = 1;
            break;
         end
      end
      if test_failed, break; end;
   end
   if test_failed, break; end;
end
function test_failed = test_multid
test_failed = 0;
h = wfilt_db(4);

maxJ = 10;
mi{1}=h{2};

for ii=2:maxJ
   mitemp = multid(h,ii); 
   for jj=1:ii-1
     err =  norm(mitemp{end+1-jj}(:)-mi{jj}(:));
     if(err>1e-6)
         test_failed = 1;
         return;
     end
   end
   mi{ii} = mitemp{2};
end    
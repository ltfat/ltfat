function test_failed=test_wfbt2filterbank
%TEST_FILTERBANK test the equivalence of wfbt anf filterbank
%  Usage: test_wfbt2filterbank()
%


disp(' ===============  TEST_WFBT2FILTERBANK ===========');

test_failed=0;

L = [ 107, 211];
W = [1,2,3];
gt = {};
gt{1} = {'db4',1,'dwt'};
gt{2} = {'db4',4,'full'};
gt{3} = {{'ana:spline4:4',3,'dwt'},{'syn:spline4:4',3,'dwt'}};
gt{4} = wfbtinit(gt{2});
gt{4} = wfbtremove(3,0,gt{4});
gt{4} = wfbtremove(3,1,gt{4});


crossoverval = 10000;

for jj=1:numel(gt)

   for ww=1:numel(W)
   for ii=1:numel(L)
      f = tester_rand(L(ii),1);
      gttmp = gt(jj);
         if iscell(gt{jj}) && iscell(gt{jj}{1})
           gttmp = gt{jj}{1};
         else
           gttmp = gttmp{1};
         end

      
      refc = wfbt(f,gttmp);

      [g,a] = wfbt2filterbank(gttmp);
      
      
      c = filterbank(f,g,a,'crossover',1);
      
      err = norm(cell2mat(c)-cell2mat(refc));
      [test_failed,fail]=ltfatdiditfail(err,test_failed);
      
      if ~isstruct(gt{jj})
         fprintf('FILT %d, COEF, FFT         L= %d, W= %d, err=%.4e %s \n',jj,L(ii),W(ww),err,fail); 
      else
         fprintf('FILT %d, COEF, FFT         L= %d, W= %d, err=%.4e %s\n',jj,L(ii),W(ww),err,fail); 
      end;
      
       c = filterbank(f,g,a,'crossover',crossoverval);
       
       err = norm(cell2mat(c)-cell2mat(refc));
       [test_failed,fail]=ltfatdiditfail(err,test_failed);
       
       if ~isstruct(gt{jj})
          fprintf('FILT %d, COEF,  TD         L= %d, W= %d, err=%.4e %s \n',jj,L(ii),W(ww),err,fail); 
       else
          fprintf('FILT %d, COEF,  TD         L= %d, W= %d, err=%.4e %s\n',jj,L(ii),W(ww),err,fail); 
       end;
      
      if iscell(gt{jj}) && iscell(gt{jj}{2})
         gttmp = gt{jj}{2};
         [g,a] = wfbt2filterbank(gttmp);
      end
      
      
      fhat = ifilterbank(refc,g,a,L(ii),'crossover',1);
      err = norm(fhat-f);
      [test_failed,fail]=ltfatdiditfail(err,test_failed);
      
      if ~isstruct(gt{jj})
         fprintf('FILT %d, INV,  FFT         L= %d, W= %d, err=%.4e %s \n',jj,L(ii),W(ww),err,fail); 
      else
         fprintf('FILT %d, INV,  FFT         L= %d, W= %d, err=%.4e %s\n',jj,L(ii),W(ww),err,fail); 
      end;

      fhat = ifilterbank(refc,g,a,L(ii),'crossover',crossoverval);
      err = norm(fhat-f);
      [test_failed,fail]=ltfatdiditfail(err,test_failed);
      
      if ~isstruct(gt{jj})
         fprintf('FILT %d, INV,   TD         L= %d, W= %d, err=%.4e %s \n',jj,L(ii),W(ww),err,fail); 
      else
         fprintf('FILT %d, INV,   TD         L= %d, W= %d, err=%.4e %s\n',jj,L(ii),W(ww),err,fail); 
      end;
      

      
      
      gttmp = gt(jj);
      if iscell(gt{jj}) && iscell(gt{jj}{1})
         gttmp = gt{jj}{1};
      else
         gttmp = gttmp{1};
      end
      
      urefc = uwfbt(f,gttmp);

      g = wfbt2filterbank(gttmp);
      
      uc = ufilterbank(f,g,1,'crossover',crossoverval);
      
      err = norm(uc-urefc);
      [test_failed,fail]=ltfatdiditfail(err,test_failed);
      
      if ~isstruct(gt{jj})
         fprintf('FILT %d, COEF, FFT, UNDEC, L= %d, W= %d, err=%.4e %s \n',jj,L(ii),W(ww),err,fail); 
      else
         fprintf('FILT %d, COEF, FFT, UNDEC, L= %d, W= %d, err=%.4e %s\n',jj,L(ii),W(ww),err,fail); 
      end;
      
      if iscell(gt{jj}) && iscell(gt{jj}{2})
         gttmp = gt{jj};
         g = wfbt2filterbank(gttmp{2});
      end
      fhat = ifilterbank(urefc,g,ones(numel(g),1),L(ii),'crossover',crossoverval);
      
      err = norm(fhat-f);
      [test_failed,fail]=ltfatdiditfail(err,test_failed);
      
      if ~isstruct(gt{jj})
         fprintf('FILT %d, INV,  FFT, UNDEC, L= %d, W= %d, err=%.4e %s \n',jj,L(ii),W(ww),err,fail); 
      else
         fprintf('FILT %d, INV,  FFT, UNDEC, L= %d, W= %d, err=%.4e %s\n',jj,L(ii),W(ww),err,fail); 
      end;
      
   end
   end
end
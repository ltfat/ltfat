function test_failed = test_wfbt(verbose)
%TEST_WFBTPR
%
% Checks perfect reconstruction of the general wavelet transform of different
% filters
%
disp('========= TEST WFBT ============');
global LTFAT_TEST_TYPE;
tolerance = 2e-8;
if strcmpi(LTFAT_TEST_TYPE,'single')
   tolerance = 1e-5;
end


test_failed = 0;
if(nargin>0)
   verbose = 1;
else
   verbose = 0;
end

type = {'dec'};
ext = {'per','zero','odd','even'};


J = 4;
%! Mild tree
 wt1 = wfbtinit({'db10',6,'full'});
 wt1 = wfbtremove(2,1,wt1,'force');
 wt1 = wfbtremove(2,3,wt1,'force');

% %! Hardcore tree
wt2 = wfbtinit({'db3',1});
wt2 = wfbtput(1,1,'mband1',wt2);
  wt2 = wfbtput(2,2,'mband1',wt2);
  wt2 = wfbtput(3,3,'mband1',wt2);
  wt2 = wfbtput(3,1,'db10',wt2);
  wt2 = wfbtput(4,1,'dgrid2',wt2);
  wt2 = wfbtput(5,1,'db3',wt2);


%! Another tree
wt3 = wfbtinit();
wt3 = wfbtput(0,0,'cmband4',wt3);
wt3 = wfbtput(1,0,'cmband6',wt3);
wt3 = wfbtput(1,1,'cmband6',wt3);
wt3 = wfbtput(2,0,'cmband4',wt3);
wt3 = wfbtput(2,1,'cmband4',wt3);
wt3 = wfbtput(3,1,'cmband4',wt3);
wt3 = wfbtput(3,2,'cmband4',wt3);
wt3 = wfbtput(3,3,'cmband4',wt3);
wt3 = wfbtput(3,4,'cmband4',wt3);

% wt2 = wfbtinit();
% wt2 = wfbtput(0,0,{'db',4},wt2);
% wt2 = wfbtput(1,0,{'algmband',1},wt2);
% wt2 = wfbtput(1,1,{'hden',3},wt2);
% wt2 = wfbtput(2,0,{'dgrid',2},wt2);
% wt2 = wfbtput(2,1,{'dgrid',2},wt2);



test_filters = {
               {'algmband2',J} % 4 filters, uniform, crit. sub.
             %  {'db4',J}
              % {'algmband1',J} % 3 filters, uniform, crit. sub.
               %{{'hden',3},J} % 3 filters, non-uniform, no crit. sub. no correct
               {'dgrid1',J} % 4 filters. sub. fac. 2
               wt1
               wt2
               wt3
               };





%testLen = 4*2^7-1;%(2^J-1);
testLen = 53;
f = tester_rand(testLen,10);

for extIdx=1:length(ext)  
   extCur = ext{extIdx};

   for typeIdx=1:length(type)
     for tt=1:length(test_filters)
        actFilt = test_filters{tt};
         if verbose, if(~isstruct(actFilt))fprintf('J=%d, filt=%s, ext=%s, inLen=%d \n',actFilt{2},actFilt{1},extCur,size(f,1)); else disp('Custom'); end; end;

        [c,info] = wfbt(f,actFilt,extCur);
        fhat = iwfbt(c,actFilt,size(f,1),extCur);
        
        %MSE
            err = norm(f-fhat,'fro');
            [test_failed,fail]=ltfatdiditfail(err,test_failed,tolerance);
            if(~verbose)
              if(~isstruct(actFilt))fprintf('J=%d, %5.5s, ext=%s, L=%d, err=%.4e %s \n',actFilt{2},actFilt{1},extCur,size(f,1),err,fail); else fprintf('Custom, err=%.4e %s\n',err,fail); end;
            end
            if strcmpi(fail,'FAILED')
               if verbose
                if(~isstruct(actFilt)) fprintf('err=%d, filt=%s, ext=%s, inLen=%d \n',err,actFilt{1},extCur,testLen); else disp('Custom'); end;
                 figure(1);clf;stem([f,fhat]);
                 figure(2);clf;stem([f-fhat]);
                 break; 
               end
            end
            
         fhat2 = iwfbt(c,info);
         err = norm(f-fhat2,'fro');
         [test_failed,fail]=ltfatdiditfail(err,test_failed,tolerance);  
         if(~isstruct(actFilt))
             fprintf('INFO J=%d, %5.5s, ext=%s, L=%d, err=%.4e %s \n',actFilt{2},actFilt{1},extCur,size(f,1),err,fail); 
         else
             fprintf('INFO Custom, err=%.4e %s\n',err,fail); 
         end;
            
            if test_failed && verbose, break; end;
        
     end
     if test_failed && verbose, break; end;
   end
   if test_failed && verbose, break; end;
end




















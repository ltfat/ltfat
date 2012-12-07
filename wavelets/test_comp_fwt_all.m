function test_failed = test_comp_fwt_all
%TEST_COMP_FWT_ALL
%
%  Checks perfect reconstruction of  comp_fwt_all comp_ifwt_all for
%  multiple channel input, usind Daubechies filters m=2,...,20, J=1,...,10
% 
%
test_failed = 0;

%  curDir = pwd;
%  mexDir = [curDir(1:strfind(pwd,'\ltfat')+5),'\mex'];
% 
%  rmpath(mexDir);
%  which comp_fwt_all
%  addpath(mexDir)
%  which comp_fwt_all
 
which comp_fwt_all -all
which comp_ifwt_all -all

type = {'dec','undec'};
ext = {'ppd','per'};
%ext = {'per','zpd','sym','symw','asym','asymw','ppd','sp0'};

f = randn(14576,1);
% multiple channels
f = [2*f,f,0.1*f];
    
for typeIdx=1:length(type)
  typeCur = type{typeIdx};
  for extIdx=1:length(ext)  
     extCur = ext{extIdx};  
     
     for ord=1:10
         [H, G, a] = wfilt_db(ord);  
         for J=1:9
            c = comp_fwt_all(f,H,J,a,typeCur,extCur);
            fhat = comp_ifwt_all(c,G,J,a,length(f),typeCur,extCur);
            err = norm(f(:)-fhat(:));
            if err>1e-6
               stem([f,fhat]); 
               test_failed = 1; 
               fprintf('J=%d, ord=%d, type=%s, ext=%s',J,ord,extCur,typeCur);
               break;
            end
         end
         if test_failed, break; end;
     end
     if test_failed, break; end;
  end 
   if test_failed, break; end;
end

 
 
   

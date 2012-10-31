function test_failed = test_comp_fwt_all
%TEST_COMP_FWT_ALL
%
% Checks perfect reconstruction of the wavelet transform
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
ext = {'per','ppd'};
%ext = {'per','zpd','sym','symw','asym','asymw','ppd','sp0'};

f = randn(14576,1);
f = [2*f,f,0.1*f];
    
for typeIdx=1:length(type)
  typeCur = type{typeIdx};
  for extIdx=1:length(ext)  
     extCur = ext{extIdx};  
     
     for ord=1:10
         for J=1:9
            [H, G] = dbfilt(ord);  
            c = comp_fwt_all(f,H,J,typeCur,extCur);
            fhat = comp_ifwt_all(c,G,J,length(f),typeCur,extCur);
            err = norm(f(:)-fhat(:));
            if err>1e-10
               test_failed = 1; 
            end
         end
     end
  end 
end

 
 
   

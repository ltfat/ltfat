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

type = {'undec'};
ext = {'ppd'};
%ext = {'per','zpd','sym','symw','asym','asymw','ppd','sp0'};

f = randn(67,1);

% multiple channels
%f = [2*f,f,0.1*f];
    
for typeIdx=1:length(type)
  typeCur = type{typeIdx};
  for extIdx=1:length(ext)  
     extCur = ext{extIdx};  
     
     for ord=1:1
         %[H, G, a] = wfilt_dden(2);  
         %[H, G, a] = wfilt_db(6); 
         %[H, G, a] = wfilt_mband(1);
         [H, G, a] = wfilt_algmband(1);
         
         for J=1:3
            c = comp_fwt_all(f,H,J,a,typeCur,extCur);
            fhat = comp_ifwt_all(c,G,J,a,length(f),typeCur,extCur);
            figure(3);clf;printCoeffs( c);
            err = norm(f(:)-fhat(:));
          % err=0;
            if err>1e-6
               figure(1);clf;stem([f,fhat]); 
               figure(2);clf;stem(f-fhat); 
               figure(3);clf;printCoeffs( c);
               test_failed = 1; 
               fprintf('error=%d,J=%d, ord=%d, type=%s, ext=%s',err,J,ord,extCur,typeCur);
               break;
            end
         end
         if test_failed, break; end;
     end
     if test_failed, break; end;
  end 
   if test_failed, break; end;
end

 

function printCoeffs( x)

[J,N1] = size(x);

for j=1:J
    subplot(J,1,j);
    stem([x{j}(:)]);
end
 
   

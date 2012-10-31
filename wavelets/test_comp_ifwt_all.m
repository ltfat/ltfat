function test_failed = test_comp_ifwt_all

% curDir = pwd;
% mexDir = [curDir(1:strfind(pwd,'\ltfat')+5),'\mex'];
% which comp_fwt_all
% rmpath(mexDir);
% which comp_fwt_all
% addpath(mexDir)
% which comp_fwt_all


type = {'dec','undec'};
ext = {'per','zpd','sym','symw','asym','asymw','ppd','sp0'};


    
for typeIdx=1:length(type)
  typeCur = type{typeIdx};
  for extIdx=1:length(ext)  
     extCur = ext{extIdx};  
     
     for ord=1:10
         for J=1:10
            [H, G] = dbfilt(ord);  
            c = comp_fwt_all(f,H,J,typeCur,extCur);
            fhat = comp_ifwt_all(c,G,J,length(f),typeCur,extCur);
            err = norm(f-fhat);
            if err>1e-10
               test_failed = 1; 
            end
         end
     end
  end 
end


   close all;
    wavName = sprintf('db%d',ord);

    dwtm = 'sp0';
    dwtmode(dwtm,'nodisp');


    %f = 0:ord*2^J+10 -1;
    %f = f(:);
    f = randn((2^J-1)*(length(H{1})-1)+5,1);
    f = [2*f];
    %f = randn(3,1);

   %tic; [C,L] = wavedec(f,J,wavName); toc;

    tic; 
        [C,L] = wavedec(f,J,H{1},H{2});
        %coefs = fwt(f,H,J);
    t1=toc;

   tic; 
      c = comp_fwt_all(f,H,J,'dec',dwtm);
   t2=toc;
   
    tic; 
        xhat1 = waverec( C,L,G{1},G{2} );
    t3=toc; 
   
   tic; 
     fhat = comp_ifwt_all(c,G,J,length(f),'dec',dwtm);
   t4=toc; 
   
   
   stem(0:length(f)-1,[f,fhat]);
   err = norm(f-fhat);
   if err>1e-10
        error('Error in reconstruction: %g',err);
   end
   
fprintf('Forward speedup: %f \n', t1/t2);
fprintf('Inverse speedup: %f \n', t3/t4);
%     if(norm(f(:)-fhat(:))>1^-10)
%      figure(2);clf;
%      stem([f,fhat]);
%     end

     [err,coefs] = checkCoefs(c,C,L,G{1},G{2},J);
%err = checkCoefsInCells(c,coefs,J);
      if(err>10^-10)
     %  figure(1);
      % clf;
     %printCoeffs(c,coefs);
           warning('Coefficients are not equal! Error is %g',err);
      end




function [err, coefs]  = checkCoefs(c,C,S,lo_r,hi_r,J)

coefs = cell(J+1,1);

coefs{1,1} = appcoef(C,S,lo_r,hi_r,J);
for j=1:J
     [coefs{end-j+1}] = detcoef(C,S,j); 
end

err = checkCoefsInCells(c,coefs,J);


function err = checkCoefsInCells(c,coefs,J)

err = 0;
err = err +  norm(c{J+1,1}(:) - coefs{J+1,1}(:));
for j=1:J
     err = err +  norm(c{j,1}(:) - coefs{j,1}(:)); 
end

function printError( x,y)

[J,N1] = size(x);

for j=1:J
    subplot(J,1,j);
     err = x{j} - y{j};
      stem(err);
      lh = line([0 length(x{j})],[eps eps]);
      set(lh,'Color',[1 0 0]);
      lh =line([0 length(x{j})],[-eps -eps]);
      set(lh,'Color',[1 0 0]);

end

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


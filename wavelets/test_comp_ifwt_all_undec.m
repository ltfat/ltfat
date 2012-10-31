function test_failed = test_comp_ifwt_all_undec

%close all;
ord = 10;
J=10; 
   close all;
    wavName = sprintf('db%d',ord);
    [H, G] = dbfilt(ord);
    dwtm = 'sym';
    dwtmode(dwtm,'nodisp');


    %f = 0:145 -1;
    %f = f(:);
    % f = randn(2^(J-1)*(length(H{1})-1),1);
    f = randn(2,1);

   %tic; [C,L] = wavedec(f,J,wavName); toc;



   tic; 
      c = comp_fwt_all(f,H,J,'undec',dwtm);
   t2=toc;
   
   
   tic; 
     fhat = comp_ifwt_all(c,G,J,length(f),'undec',dwtm);
   t4=toc; 
   
    err = norm(f-fhat);
 stem(0:length(f)-1,[f, fhat]);

   if err>10^(-10)
        error('Error in reconstruction: %g',err);
   end

   
%fprintf('Forward speedup: %f \n', t1/t2);
%fprintf('Inverse speedup: %f \n', t3/t4);
%     if(norm(f(:)-fhat(:))>1^-10)
%      figure(2);clf;
%      stem([f,fhat]);
%     end

%      [err,coefs] = checkCoefs(c,C,L,G{1},G{2},J);
% %err = checkCoefsInCells(c,coefs,J);
%       if(err>1^-10)
%      %  figure(1);
%        clf;
%      %  printCoeffs(c,coefs);
%            error('Coefficients are not equal! Error is %g',err);
%       end




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


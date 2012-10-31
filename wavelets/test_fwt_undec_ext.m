function test_failed = test_fwt_undec_ext

close all;
    ord = 2;
    J= 4;
   
    wavName = sprintf('db%d',ord);
    [H, G] = dbfilt(ord);
    dwtM = 'ppd';
    dwtmode(dwtM,'nodisp');


    f = 1:2^J*10;
    f = f(:);
  %  f = greasy;
   % f = randn(2^J*180+47,1);
   % f= zeros(2^J*180,1);
f = [randn(2^J*10,1)];
    

  %tic; 
    c = comp_fwt_all(f,H,J,'undec',dwtM);
  %t2=toc;
 % fprintf('Speedup: %f \n', t1/t2);
    %fhat = ifwt(c,G,J,length(f),'undec');
    
%     if(norm(f-fhat)>1^-10)
%     figure(2); clf;
%     stem([f,fhat]);
%     end
  %  coefs = W.dec.';
 %    err = checkCoefsInCells(c,coefs,J);
      printCoeffs(c,c);
     if(err>1^-10)
         figure(1);clf;
         printCoeffs(c,coefs);
         error('Coefficients are not equal! Error is %g',err);
    end




function [err, coefs]  = checkCoefs(c,SWC,J)

coefs = cell(J+1,1);

coefs{1,1} = SWC(end,:);
for j=1:J
     [coefs{end-j+1}] = SWC(j,:); 
end


err = 0;
err = err +  norm(c{J+1,1}(:) - coefs{J+1,1}(:));
for j=1:J
     err = err +  norm(c{j,1}(:) - coefs{j,1}(:)); 
end

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
     err = x{j}(:) - y{j}(:);
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




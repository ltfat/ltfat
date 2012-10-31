function hat2 = test_comp_up_conv_td

%close all;
ord = 8;
J=4; 
   close all;
    wavName = sprintf('db%d',ord);
    [H, G] = dbfilt(ord);
    fLen = length(G{1});
    skip = 0;
    dwtm = 'ppd';
    dwtmode(dwtm,'nodisp');
    filtUps = 0;
    up = 0;
    inLen = 101;
    %f = 0:inLen -1;
    %f = {f(:);f(:) };

    f ={ randn(inLen,1);  randn(inLen,1)};
    outLen =inLen+1000;%up*inLen+fLen;
    
   upsf1 = ups(f{1},up+1,3);
   upsf2 = ups(f{2},up+1,3);
   fhat1 = conv(upsf1,ups(G{1},2^filtUps,1)) + conv(upsf2,ups(G{2},2^filtUps,1));
   skip = (2^filtUps)*fLen - (2^filtUps-1) - 1;
   fhat1 = fhat1(1+skip:min([outLen,length(fhat1) ]));
   
   fhat2 = comp_up_conv_td(f,G,length(fhat1),up,skip,filtUps);
   
   stem(0:length(fhat1)-1,[fhat1,fhat2]);
   err = norm(fhat1-fhat2);
   if err>1^-10
        error('Error in reconstruction: %g',err);
   end
   
% fprintf('Speedup: %f \n', t1/t2);
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


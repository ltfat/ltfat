function test_failed = test_conv_sub


ord = 6;
J=1;
   
    wavName = sprintf('db%d',ord);
    [H, G] = dbfilt(ord);
    dwtmode('zpd','nodisp');


    %f = 0:157;
    %f = f(:);
    f = randn(175,1);

   %tic; [C,L] = wavedec(f,J,wavName); toc;
sub = 2;
filtSub = 1;
coefs= cell(J+1,1);
 tic;
  coefs{1} = downs(conv(f,H{1}),sub);
  coefs{2} = downs(conv(f,H{2}),sub);
 t1=toc;


htemp = cell(1);
htemp{1} = H{1};
   tic;  c = comp_fwt_all(f,H,sub,sub-1,0,filtSub); t2=toc;
fprintf('Speedup: %f \n', t1/t2);
%     if(norm(f(:)-fhat(:))>1^-10)
%      figure(2);clf;
%      stem([f,fhat]);
%     end

 %    [err,coefs] = checkCoefs(c,C,L,G{1},G{2},J);
err = checkCoefsInCells(c,coefs,J);
      if(err>1^-10)
       figure(1);
       clf;
       printCoeffs(c,coefs);
           error('Coefficients are not equal! Error is %g',err);
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


function test_failed = test_fwt_vs_wtfft(verbose)
%TEST_COMP_FWT_ALL
%
% Checks equality of coefficients obtained by functions fwt and wtfft and
% perfect reconstruction using both ifwt and iwtfft
%
test_failed = 0;
if(nargin>0)
   verbose = 1;
else
   verbose = 0;
end

which comp_fwt_all -all
which comp_ifwt_all -all


test_filters = {
               %{'lemaire',80} % not exact reconstruction
               %{'hden',3} % only one with 3 filters and different sub
               %facts
               %subsampling factors, not working
               %{'algmband',2} % 4 filters, subsampling by the factor of 4!, really long identical lowpass filter
               {'sym',4}
               {'sym',9}
               {'symds',1}
               {'symds',2}
               {'symds',3}
               {'symds',4}
               {'symds',5}
               {'spline',4,4}
               {'spline',3,5}
               {'spline',3,11}
               {'spline',11,3}
               {'maxflat',2}
               {'maxflat',11}
               {'db',1}
               {'db',10}
               {'db',20}
               %{'optfs',7} % bad precision of filter coefficients
               %{'remez',20,10,0.1} % no perfect reconstruction
               {'dden',1}
               {'dden',2}
               {'dden',3}
               {'dden',4}
               {'dden',5}
               {'dden',6}
               {'dden',7}
               {'dgrid',1}
               {'dgrid',2}
               {'dgrid',3}
               {'algmband',2} 
               {'mband',1}
               %{'hden',3}
               %{'hden',2}
               %{'hden',1}
               %{'apr',1} % complex filter values
               };

           
J = 3;
testLen = 200*2^J; % can be compared only for lengths equalt to the multiple of power of two 
f = randn(testLen,1);

for tt=1:length(test_filters)
         actFilt = test_filters{tt};
         %fname = strcat(prefix,actFilt{1});
         w = waveletfb(actFilt); 

     for jj=1:J
           if verbose, fprintf('J=%d, filt=%s \n',jj,actFilt{1}); end; 
           
             [h,a] = multid(w.h,jj,w.a);
             [g,a] = multid(w.g,jj,w.a,'syn');
             H = freqzfb(h,filterbanklength(length(f),a));
             G = freqzfb(g,filterbanklength(length(f),a));
         
             c1 = wtfft(f,H,a);
             c2 = fwt(f,w,jj);

             fhat = iwtfft(c1,G,a,length(f));
             fhat2 = ifwt(c2,w,jj,length(f));
           
             % check reconstruction
             err = norm(f(:)-fhat(:));
             err = err + norm(f(:)-fhat2(:));
            if err>1e-6
               test_failed = 1; 
               if verbose
                 fprintf('Bad reconstruction: err=%d, J=%d, filt=%s',err,jj,actFilt{1},extCur,typeCur);
                 figure(1); stem([f,fhat,fhat2]);
               end
               break; 
            end
            
            % check coefficients
            err = checkCoefsInCells(c1,changeFormat(c2,jj,length(w.h)),jj);
            if err>1e-6
               test_failed = 1; 
               if verbose
                 fprintf('Not Equal Coefficients: err=%d, J=%d, filt=%s \n',err,jj,actFilt{1});
                 printCoeffs( c1,changeFormat(c2,jj,length(w.h)));
               end
               break; 
            end
            
     end
      if test_failed, break; end;
end 






function err = checkCoefsInCells(c,coefs,J)

err = 0;
err = err +  norm(c{J+1}(:) - coefs{J+1}(:));
for j=1:J
     err = err +  norm(c{j}(:) - coefs{j}(:)); 
end


function c = changeFormat(c2,J,noFilt)
c2form = cell(numel(c2)-(noFilt-1),1);
c2form{1} = c2{1};
cSformIdx = 2;
for jj=2:J+1
    for ii=1:noFilt-1
       c2form{cSformIdx} = c2{jj,ii};
       cSformIdx=cSformIdx+1;
    end
end
c = c2form;


 
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

function coefs = coefMatToLTFAT(C,S,lo_r,hi_r,J)

coefs = cell(J+1,1);

coefs{1,1} = appcoef(C,S,lo_r,hi_r,J);
for j=1:J
     [coefs{end-j+1}] = detcoef(C,S,j); 
end


 
   

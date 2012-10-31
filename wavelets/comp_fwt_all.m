function c = comp_fwt_all(f,h,J,type,ext)
%COMP_FWT_ALL Compute DWT
%   Usage:  c=comp_fwt_all(f,h,J,type,ext);
%
%   Input parameters:
%         f     : Input data.
%         h     : Analysis Wavelet filters.
%         J     : Number of filterbank iterations.
%         type  : 'dec','undec' Type of the wavelet transform.
%         ext   : 'per','zpd','sym','symw','asym','asymw','ppd','sp0' Type of the forward transform boundary handling.
%
%   Output parameters:
%         c     : Coefficients stored in J+1 cell-array.
%

% Do non-expansve transform if ext='per'
if(strcmp(ext,'per'))
    doNoExt = 1;
else
    doNoExt = 0;
end

c = cell(J+1,1);
cLen = zeros(J,1);
fLen = length(h{1});
[inLen, chans] = size(f);


if(strcmp(type,'dec'))
   sub = 2;
   if(doNoExt)
       ext='perdec';
       skip = ceil(fLen/2);

        for jj=1:J
           cLen(J+1-jj) = ceil(2^(-jj)*inLen);
        end
   else
       skip = 1;
       for jj=1:J
           cLen(J+1-jj) = floor(2^(-jj)*inLen + (1-2^(-jj))*(fLen-1));
       end
      
   end
   
   for ch=1:chans
    tempca = f(:,ch);
      for jj=1:J
        ctemp = conv_td_sub(tempca,cLen(J+1-jj),h,sub,skip,ext,0);
        tempca = ctemp{1};
        c{J+2-jj}(:,ch) = ctemp{2};
      end
     c{1}(:,ch) = tempca;
   end
    
elseif(strcmp(type,'undec'))
   sub = 1;
   h{1} = h{1}/sqrt(2);
   h{2} = h{2}/sqrt(2);
   
   if(doNoExt)
       
    for ch=1:chans
    tempca = f(:,ch);  
      for jj=1:J
        skip = floor((2^(jj-1)*fLen+1)/2);  
        ctemp = conv_td_sub(tempca,inLen,h,sub,skip,ext,2^(jj-1));
        tempca = ctemp{1};
        c{J+2-jj}(:,ch) = ctemp{2};
      end
     c{1}(:,ch) = tempca;
    end   
        
   else
       cLen(J) = inLen + fLen -1;
       for jj=J-1:-1:1
           cLen(jj) = cLen(jj+1) + 2^(J-jj)*fLen-(2^(J-jj)-1) -1;
       end
       
    for ch=1:chans
    tempca = f(:,ch);  
      for jj=1:J
        ctemp = conv_td_sub(tempca,cLen(J+1-jj),h,sub,0,ext,2^(jj-1));
        tempca = ctemp{1};
        c{J+2-jj}(:,ch) = ctemp{2};
      end
     c{1}(:,ch) = tempca;
    end  
   end   
end


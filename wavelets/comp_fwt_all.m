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

[filts, filtCols] = size(h);
c = cell(J+1,filts-1);

cLen = zeros(J,1);
fLen = length(h{1});
[inLen, chans] = size(f);

%
% The decimated case
if(strcmp(type,'dec') || strcmp(type,'dtdwt') || strcmp(type,'hddwt'))
   sub = 2;
   if(doNoExt)
       ext='perdec';
       skip = floor(fLen/2);

        for jj=1:J
           cLen(J+1-jj) = ceil(2^(-jj)*inLen);
        end
   else
       skip = 1;
       for jj=1:J
           cLen(J+1-jj) = floor(2^(-jj)*inLen + (1-2^(-jj))*(fLen-1));
       end
   end
   
   if(strcmp(type,'dec'))
       for ch=1:chans
        tempca = f(:,ch);
          for jj=1:J
            ctemp = conv_td_sub(tempca,cLen(J+1-jj),h,sub,skip,ext,0);
            tempca = ctemp{1};
            for ff=1:filts-1
              c{J+2-jj,ff}(:,ch) = ctemp{1+ff};
            end
          end
         c{1}(:,ch) = tempca;
       end
   elseif(strcmp(type,'dtdwt'))
      for ii = 1:4
         h{ii} = h{ii}/sqrt(2);
      end
        % for all input channels
        for ch=1:chans
            tempca = f(:,ch);
            % first level, different filters
            ctemp1 = conv_td_sub(tempca,cLen(J),{h{:,1}},sub,skip,ext,0);
            ctemp2 = conv_td_sub(tempca,cLen(J),{h{:,2}},sub,skip,ext,0);
            tempca1 = ctemp1{1};
            tempca2 = ctemp2{1};
            c{J+1}(:,ch) = ctemp1{2}+i.*ctemp2{2};
              for jj=2:J
                ctemp1 = conv_td_sub(tempca1,cLen(J+1-jj),{h{:,3+mod(jj,2)}},sub,skip,ext,0);
                ctemp2 = conv_td_sub(tempca2,cLen(J+1-jj),{h{:,4-mod(jj,2)}},sub,skip,ext,0);
                tempca1 = ctemp1{1};
                tempca2 = ctemp2{1};
                for ff=1:filts-1
                  c{J+2-jj,ff}(:,ch) = ctemp1{1+ff}+i.*ctemp2{1+ff};
                end
              end
             c{1}(:,ch) = tempca1 +i*tempca2;
        end
    elseif(strcmp(type,'hddwt'))
      ch3cLen = [cLen; inLen];  
      for ch=1:chans
        tempca = f(:,ch);
          for jj=1:J
            ctemp = conv_td_sub(tempca,cLen(J+1-jj),{h{1},h{2}},sub,skip,ext,0);
            c{J+2-jj,2}(:,ch) = conv_td_sub(tempca,ch3cLen(end+1-jj),{h{3}},1,skip,ext,0);
            c{J+2-jj,1}(:,ch) = ctemp{2};
            tempca = ctemp{1};
          end
         c{1}(:,ch) = tempca;
       end
    end
% time-invariant wavelet transform    
elseif(strcmp(type,'undec'))
   sub = 1;
   for ii = 1:numel(h)
     h{ii} = h{ii}/sqrt(2);
   end
   
   if(doNoExt)
       
    for ch=1:chans
    tempca = f(:,ch);  
      for jj=1:J
        skip = ceil((2^(jj-1)*fLen)/2);  
        ctemp = conv_td_sub(tempca,inLen,h,sub,skip,ext,2^(jj-1));
        tempca = ctemp{1};
        for ff=1:filts-1
          c{J+2-jj,ff}(:,ch) = ctemp{1+ff};
        end
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
        for ff=1:filts-1
          c{J+2-jj,ff}(:,ch) = ctemp{1+ff};
        end
      end
     c{1}(:,ch) = tempca;
    end  
   end   
end


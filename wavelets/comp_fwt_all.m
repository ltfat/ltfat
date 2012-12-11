function c = comp_fwt_all(f,h,J,a,type,ext)
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

[filts] = length(h);
c = cell((filts-1)*J+1,1);

% number of coefficients at each level, according to the filtertree lowpass path
cLen = zeros(J,1);
% assuming all filters are equal length
fLen = length(h{1});
[inLen, chans] = size(f);

%
% The decimated case
if(strcmp(type,'dec'))
sub = a(1);    
if(doNoExt)
   skip = floor(fLen/2); 
   ext='perdec';
   for jj=1:J
       cLen(J+1-jj) = ceil(sub^(-jj)*inLen);
   end
else
   skip = 1; 
   cLen(J) = ceil((inLen + fLen-1 -skip)/sub);
   for jj=2:J
      cLen(J+1-jj) = ceil((cLen(J+1-jj+1)+ fLen-1-skip)/sub);%floor(sub^(-jj)*inLen + (1-sub^(-jj))*(fLen-1));
    end
end
   % all subsampling factors are equal 
   if all(a == a(1))
         for ch=1:chans
            tempca = f(:,ch);
               for jj=1:J
                  ctemp = conv_td_sub(tempca,cLen(J+1-jj),h,sub,skip,ext,0);
                  tempca = ctemp{1};
                  for ff=1:filts-1
                     c{end-jj*(filts-1)+ff}(:,ch) = ctemp{1+ff};
                  end
               end
            c{1}(:,ch) = tempca;
         end
   else
      % case in which not all subsampling factors are equal, filters are
      % taken individually
      % array of lengths of input signals
      cLen = [cLen;inLen];
      for ch=1:chans
            tempca = f(:,ch);
               for jj=1:J
                  actInLen = cLen(end+1-jj); 
                  for ff=1:filts-1
                     if(doNoExt), actOutLen = ceil(actInLen/a(ff+1));
                     else actOutLen = floor((actInLen+(fLen-1))/a(ff+1)); end 
                     
                     c{end-jj*(filts-1)+ff}(:,ch) = conv_td_sub(tempca,actOutLen,{h{ff+1}},a(ff+1),skip,ext,0);
                  end
                     tempca = conv_td_sub(tempca,cLen(J+1-jj),{h{1}},a(1),skip,ext,0);
               end
           c{1}(:,ch) = tempca;
      end
   end

% time-invariant wavelet transform    
elseif(strcmp(type,'undec'))
   sub = 1;
   % Since no downsampling takes place, normalize impulse responses
   for ii = 1:numel(h)
     h{ii} = h{ii}/sqrt(a(ii));
   end
   
    skip = zeros(J,1);
    cLen = zeros(J,1);
    if(doNoExt)
       cLen(:) = inLen;
       for jj=1:J
           skip(jj) = ceil((a(1)^(jj-1)*fLen)/2);
       end
    else
       skip(:) = 0; 
       cLen(J) = inLen + fLen -1;
       for jj=J-1:-1:1
           cLen(jj) = cLen(jj+1) + a(1)^(J-jj)*fLen-(a(1)^(J-jj)-1) -1;
       end
    end
    
    for ch=1:chans
    tempca = f(:,ch);  
      for jj=1:J
        ctemp = conv_td_sub(tempca,cLen(J+1-jj),h,sub,skip(jj),ext,a(1)^(jj-1));
        tempca = ctemp{1};
        %c{J+2-jj}(:,ch) = ctemp{2};
        for ff=1:filts-1
          c{end-jj*(filts-1)+ff}(:,ch) = ctemp{1+ff};
        end
      end
     c{1}(:,ch) = tempca;
    end 
   
%    if(doNoExt)
%     for ch=1:chans
%     tempca = f(:,ch);  
%       for jj=1:J
%         skip = ceil((2^(jj-1)*fLen)/2);  
%         ctemp = conv_td_sub(tempca,inLen,h,sub,skip,ext,2^(jj-1));
%         tempca = ctemp{1};
%         for ff=1:filts-1
%           c{J+2-jj,ff}(:,ch) = ctemp{1+ff};
%         end
%       end
%      c{1}(:,ch) = tempca;
%     end   
%         
%    else
%        cLen(J) = inLen + fLen -1;
%        for jj=J-1:-1:1
%            cLen(jj) = cLen(jj+1) + 2^(J-jj)*fLen-(2^(J-jj)-1) -1;
%        end
       
    
   %end   
end


function c=comp_fwt_undec(f,h,J)
%COMP_FWT_UNDEC  Computes undecimated Fast Wavelet Transform 
%   Usage:  c=comp_fwt_undec(f,h,J);
%  
%   Input parameters:
%         f     : Input data.
%         h     : Wavelet filters.
%         J     : Number of filterbank iterations.
%   Output parameters:
%         c     : Coefficients stored in J+1 cell-array.
%
%   See also:
%
%   Demos:
%
%   References:

%   AUTHOR : Zdenek Prusa.
%   TESTING: COMP_FWT
%   REFERENCE: COMP_FWT

[L,W] = size(f);
c = cell(J+1,1);


flen = length(h{1});
upfLen = @(j) max([flen, 2^j*flen-(2^j-1)]);


% prepare impulse responses
    for hidx=1:length(h)
        for j=0:J-1
                flenTemp = upfLen(j);
                hpad{j+1,hidx}=zeros(L,1);
                range = 1:flenTemp;
                hpad{j+1,hidx}(range) = ups(h{hidx}(:),2^j,1);
 
                % circshift to compensate for the delay
                hpad{j+1,hidx} = circshift(hpad{j+1,hidx},-ceil((2^j*flen)/2) );
         end
     end
 
     
     for w=1:W    
         tempa = f(:,w);
         for j=1:J
             c{end+1-j}(:,w) = pconv(tempa,hpad{j,2});
             tempa = pconv(tempa,hpad{j,1});
         end
         c{1}(:,w) = tempa;
     end 
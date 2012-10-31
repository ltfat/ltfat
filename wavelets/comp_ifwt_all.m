function f = comp_ifwt_all(c,g,J,Ls,type,ext)
%COMP_IFWT_ALL Compute Inverse DWT
%   Usage:  f = comp_ifwt_all(c,g,J,Ls,type,ext);
%
%   Input parameters:
%         c     : Coefficients stored in J+1 cell-array.
%         g     : Synthesis wavelet filters.
%         Ls    : Length of the reconstructed signal.
%         type  : 'dec','undec' Type of the wavelet transform.
%         ext   : 'per','zpd','sym','symw','asym','asymw','ppd','sp0' Type of the forward transform boundary handling.
%
%   Output parameters:
%         f     : Reconstructed data.
%


if(strcmp(ext,'per'))
    doNoExt = 1;
else
    doNoExt = 0;
end

fLen = length(g{1});
len = zeros(J+1,1);
    for jj=1:J
       len(jj) = length(c{jj+1});
    end
    len(jj+1) = Ls;

[cLen,chans] = size(c{end});   
f = zeros(Ls,chans);


if(strcmp(type,'dec'))
    upFac= 2;
    if(doNoExt)
        skip = floor((fLen+1)/2) -1;
    else
        skip = fLen - 2;
    end

    for ch=1:chans
       tempca = c{1}(:,ch);
       for jj=2:J+1
           tempca = up_conv_td({tempca;c{jj}(:,ch)},len(jj),g,upFac,skip,doNoExt,0);
       end
       f(:,ch) = tempca;
    end
    

elseif(strcmp(type,'undec'))
    upFac= 1;
    g{1} = g{1}/sqrt(2);
    g{2} = g{2}/sqrt(2);
    
    if(doNoExt)
       for ch=1:chans
          tempca = c{1}(:,ch);
          for jj=2:J+1
             filtUps = 2^(J+1-jj); 
             skip = floor((filtUps*fLen)/2) - filtUps;  
             tempca = up_conv_td({tempca;c{jj}(:,ch)},len(jj),g,upFac,skip,doNoExt,filtUps);
          end
          f(:,ch) = tempca;
       end
    else
       for ch=1:chans
          tempca = c{1}(:,ch);
          for jj=2:J+1
             filtUps = 2^(J+1-jj); 
             skip = filtUps*fLen - (filtUps-1) - 1;  
             tempca = up_conv_td({tempca;c{jj}(:,ch)},len(jj),g,upFac,skip,doNoExt,filtUps);
          end
          f(:,ch) = tempca;
       end
    end
end
    
    
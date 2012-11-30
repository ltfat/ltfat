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

filts = numel(g);
fLen = length(g{1});
len = zeros(J+1,1);
    for jj=1:J
       len(jj) = length(c{jj+1});
    end
    len(jj+1) = Ls;

[cLen,chans] = size(c{end});   
f = zeros(Ls,chans);
tmpin = cell(filts,1);

if(strcmp(type,'dec') || strcmp(type,'dtdwt') || strcmp(type,'hddwt'))
    upFac= 2;
    if(doNoExt)
        skip = floor((fLen)/2) -1;
    else
        skip = fLen - 2;
    end

 if(strcmp(type,'dec'))
    for ch=1:chans
       tempca = c{1}(:,ch);
       for jj=2:J+1
           tmpin{1} = tempca;
           for ff=1:filts-1
               tmpin{1+ff}= c{jj,ff}(:,ch);
           end
           tempca = up_conv_td(tmpin, len(jj),g,upFac,skip,doNoExt,0);
       end
       f(:,ch) = tempca;
    end
 elseif(strcmp(type,'dtdwt'))
    for ii = 1:4
       g{ii} = g{ii}/sqrt(2);
    end
    for ch=1:chans
       tempca1 = real(c{1}(:,ch));
       tempca2 = imag(c{1}(:,ch));
       
       for jj=2:J
           tmpin{1} = tempca1;
           tmpin{2} = real(c{jj}(:,ch));
           tempca1 = up_conv_td(tmpin, len(jj),{g{:,3+mod(J+2-jj,2)}},upFac,skip,doNoExt,0);
           tmpin{1} = tempca2;
           tmpin{2} = imag(c{jj}(:,ch));
           tempca2 = up_conv_td(tmpin, len(jj),{g{:,4-mod(J+2-jj,2)}},upFac,skip,doNoExt,0);
       end
       tmpin{1} = tempca1;
       tmpin{2} = real(c{end}(:,ch));
       tempca1 = up_conv_td(tmpin, len(end),{g{:,1}},upFac,skip,doNoExt,0);
       tmpin{1} = tempca2;
       tmpin{2} = imag(c{end}(:,ch));
       tempca2 = up_conv_td(tmpin, len(end),{g{:,2}},upFac,skip,doNoExt,0);
       
       f(:,ch) = (tempca1 + tempca2);
    end 
 elseif(strcmp(type,'hddwt'))
     skip = floor((fLen)/2)-1;
     %skipch3 = floor((fLen+1)/2);
     for ch=1:chans
       tempca = c{1}(:,ch);
       for jj=2:J+1
           tmpin{1} = tempca;
           tmpin{2}= c{jj,1}(:,ch);
           tempca = up_conv_td(tmpin, len(jj),{g{1},g{2}},upFac,skip,doNoExt,0);
           tempca = tempca + up_conv_td({c{jj,2}(:,ch)}, len(jj),{g{3}},1,skip,doNoExt,0);
       end
       f(:,ch) = tempca;
    end
 end

elseif(strcmp(type,'undec'))
    upFac= 1;
    for ii = 1:numel(g)
       g{ii} = g{ii}/sqrt(2);
    end
    
    if(doNoExt)

       for ch=1:chans
          tempca = c{1}(:,ch);
          for jj=2:J+1
             tmpin{1} = tempca;
             for ff=1:filts-1
                 tmpin{1+ff}= c{jj,ff}(:,ch);
             end
             filtUps = 2^(J+1-jj); 
             skip = floor((filtUps*fLen)/2) - filtUps;  
             tempca = up_conv_td(tmpin,len(jj),g,upFac,skip,doNoExt,filtUps);
          end
          f(:,ch) = tempca;
       end
    else
       for ch=1:chans
          tempca = c{1}(:,ch);
          for jj=2:J+1
             tmpin{1} = tempca;
             for ff=1:filts-1
               tmpin{1+ff}= c{jj,ff}(:,ch);
             end
             filtUps = 2^(J+1-jj); 
             skip = filtUps*fLen - (filtUps-1) - 1;  
             tempca = up_conv_td(tmpin,len(jj),g,upFac,skip,doNoExt,filtUps);
          end
          f(:,ch) = tempca;
       end
    end
end
    
    
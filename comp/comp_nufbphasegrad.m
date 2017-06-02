function [tgrad,fgrad,logs] = comp_nufbphasegrad(abss,N,a,M,sqtfr,fc,NEIGH,posInfo,gderivweight)
%chanStart = [0;cumsum(N)];
fac = gderivweight;
%cfreqdiff = diff(fc);
%sqtfr = sqrt(tfr);
%sqtfrdiff = diff(sqtfr);

L = a(1)*N(1);

difforder = 2;
tt=-11;

logs = log(abss + realmin);
%logsMax = max(logs);
%logs(logs<logsMax+tt) = tt;

% Obtain the (relative) phase difference in frequency direction by taking
% the time derivative of the log magnitude and weighting it by the
% time-frequency ratio of the appropriate filter.
% ! Note: This disregards the 'quadratic' factor in the equation for the 
% phase derivative !

%tmagdiff = zeros(size(logs));
fgrad = zeros(size(logs));
chanStart = 0;
for kk = 1:M
    idx = chanStart+(1:N(kk));
    fgrad(idx) = pderiv(logs(idx),1,difforder)/N(kk);
   % fgrad(idx) = tmagdiff(idx).*sqtfr(kk)^2/(2*pi);
    
    %tmagdiff(idx) = tmagdiff(idx)/N(kk);
    
    chanStart = chanStart + N(kk);
end

% Obtain the (relative) phase difference in time direction using the
% frequency derivative of the log magnitude. The result is the mean of
% estimates obtained from 'above' and 'below', appropriately weighted by
% the channel distance and the inverse time-frequency ratio of the
% appropriate filter.
% ! Note: We consider the term depending on the time-frequency ratio 
% difference, but again disregard the 'quadratic' factor. !
%fac = 0;
%fac = 1/2; 
%fac = 2/3;
%fac = 2/pi;

tgrad = zeros(size(abss));

chanStart = 0;
for kk = 1:M
    aboveNom = 0; aboveDenom = 1; belowNom = 0; belowDenom = 1; 
    denom = sqtfr(kk)^2*(pi*L);
    if kk<M
        aboveNom = fac*(sqtfr(kk+1)-sqtfr(kk))/sqtfr(kk);
        aboveDenom = fc(kk+1)-fc(kk);
    end
    if kk>1
        belowNom = fac*(sqtfr(kk)-sqtfr(kk-1))/sqtfr(kk);
        belowDenom = fc(kk)-fc(kk-1);
    end
   
    temp = zeros(N(kk),1);    
    for ll = 1:N(kk) 
        w = chanStart + ll;
        tempValAbove = 0;
        numNeigh = 0;
        for jj = 1:2
           neigh = NEIGH(4+jj,w);           
           if neigh
              numNeigh = numNeigh+1;
              dist = (posInfo(neigh,2)-posInfo(w,2))/a(kk);
              tempValAbove = tempValAbove + (logs(neigh)-logs(w) - dist*fgrad(w));
           end
        end
        if numNeigh
           tempValAbove = tempValAbove/numNeigh;
        end
        
        tempValBelow = 0;
        numNeigh = 0;
        for jj = 1:2
           neigh = NEIGH(2+jj,w);           
           if neigh
              numNeigh = numNeigh+1;
              dist = (posInfo(neigh,2)-posInfo(w,2))/a(kk);
              tempValBelow = tempValBelow + (logs(w)-logs(neigh) - dist*fgrad(w));
           end
        end
        
        if numNeigh
            tempValBelow = tempValBelow/numNeigh; 
        end  
        
        temp(ll) = (tempValAbove + aboveNom) / aboveDenom + ...
                   (tempValBelow + belowNom) / belowDenom;
        %temp(ll,2) = (tempValBelow + belowNom) / belowDenom;
    end
%     if kk<M
%         temp(:,1) = (temp(:,1) + fac*(sqtfr(kk+1)-sqtfr(kk))./sqtfr(kk))./(fc(kk+1)-fc(kk));
%     end
%     if kk>1
%         temp(:,2) = (temp(:,2) + fac*(sqtfr(kk)-sqtfr(kk-1))./sqtfr(kk))./(fc(kk)-fc(kk-1));
%     end
    % Maybe a factor of 1/2 is missing here?
    
    tgrad(chanStart+(1:N(kk))) = temp/denom;
    
    chanStart = chanStart + N(kk);
end

chanStart = 0;
for kk = 1:M
    idx = chanStart+(1:N(kk));
    fgrad(idx) = fgrad(idx).*sqtfr(kk)^2/(2*pi)*N(kk);
    
    chanStart = chanStart + N(kk);
end




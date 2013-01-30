function H = wtfftfreqz(h,varargin);
%FREQZFB Zero-delay frequency responses from the impulse responses
%   Usage:  freqzfb(h,
%
%   XXX Merge this function with magresp
%


definput.keyvals.Ls = [];
[flags,kv,Ls]=ltfatarghelper({'Ls'},definput,varargin);   

if(isempty(Ls))
   % find the longest filter
   maxLen = length(h{1}.h);
   for ff = 2:length(h)
      tmpLen = length(h{ff}.h);
      if(tmpLen>maxLen)
         maxLen = tmpLen;
      end
   end 
   % if bigger than 1024, use next fast fft length
   Ls = nextfastfft(max([maxLen,1024])); 
end

if(isnumeric(h))
    h={h};
end

fNo = numel(h);


H = cell(fNo,1);
for jj=1:fNo
   hlen = length(h{jj}.h);
   H{jj} = wfiltstruct('FREQ');
   H{jj}.d = 0;
   if(hlen<Ls) 
      %well behaved case
      htemp = zeros(Ls,1);
      htemp(1:hlen) = h{jj}.h(:);
      H{jj}.H = fft(circshift(htemp,-(h{jj}.d-1)));
   else
      %periodic wraparound od the filter impulse response
      ratio = hlen/Ls;
      htemp = zeros(1,ceil(ratio)*Ls);
      htemp(1:hlen) = h{jj}.h(:);
      htemp = circshift(htemp,-(h{jj}.d-1));
      H{jj}.H = fft(sum(reshape(htemp,ceil(ratio),Ls)));
   end
end 


% DO THE PLOTTING
if(nargout==0)
    Hplot = zeros(Ls,fNo);
    for ff=1:fNo
       Hplot(:,ff) = H{ff}.H; 
    end
    
    
    f = -ceil((Ls-1)/2):floor((Ls-1)/2);
    f = f/(Ls/2);
    Hplot = fftshift(Hplot);
    
    % plot phase frequency response
    subplot(2,1,2);
    plot(f,unwrap(angle(Hplot))/pi);
    axis tight;  
    
    flags.do_norm = 1;
    flags.do_db = 0;

    if(flags.do_db)
       Hplot = 20*log10(abs(Hplot));
    else
       Hplot = abs(Hplot); 
    end
    
    if(flags.do_norm)
       for nn = 1:fNo
          Hplot(:,nn) = Hplot(:,nn)/max(abs(Hplot(:,nn)));
       end
    end
    
    
    
    % plot 
    subplot(2,1,1);
    plot(f,Hplot);
    axis tight;

end;
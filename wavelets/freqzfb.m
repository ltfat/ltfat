function [H,f] = freqzfb(h,varargin);
%FREQZFB Zero-delay frequency responses from the impulse responses
%   Usage:  freqzfb(h,
%
%   XXX Merge this function with magresp


definput.keyvals.Ls = 4*1024;
[flags,kv,Ls]=ltfatarghelper({'Ls'},definput,varargin);   

if(isnumeric(h))
    h={h};
end

fNo = numel(h);
H = zeros(Ls,fNo);

htemp = zeros(Ls,1); 

  for jj=1:fNo
      htemp(:) = 0;
      htemp(1:length(h{jj})) = h{jj}(:);
      H(:,jj) = fft(circshift(htemp,-floor(length(h{jj})/2)));
  end 



% DO THE PLOTTING
if(nargout==0)
    f = -ceil((Ls-1)/2):floor((Ls-1)/2);
    f = f/(Ls/2);
    H = fftshift(H);
    flags.do_norm = 1;
    flags.do_db = 0;

    if(flags.do_db)
       Hplot = 20*log10(abs(H));
    else
       Hplot = abs(H); 
    end
    
    
    
    if(flags.do_norm)
       for nn = 1:fNo
          Hplot(:,nn) = Hplot(:,nn)/max(abs(Hplot(:,nn)));
       end
    end
    
    
    

    subplot(2,1,1);
    plot(f,Hplot);
    axis tight;
    subplot(2,1,2);
    plot(f,unwrap(angle(H))/pi);
    axis tight;
end;
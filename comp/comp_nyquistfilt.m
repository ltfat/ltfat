function H = comp_nyquistfilt(wintype,fs,chan_max,freqtoscale,scaletofreq,bwmul,bins,Ls)
%COMP_NYQUISTFILT high-pass filter for warped filter banks

    kk = chan_max;
    while scaletofreq(kk-bwmul) < fs/2;
      kk = kk+1/bins;
    end
    Maxfilt = kk;
    
    Minpos = ceil(Ls/fs*scaletofreq(chan_max+1/bins-bwmul));
    samples = freqtoscale((Minpos-1:floor(Ls/2))*fs/Ls);
    
    FILTS = zeros(round(bins*(Maxfilt-chan_max)),numel(samples));
    for kk = 1:size(FILTS,1)
       FILTS(kk,:) = firwin(wintype,(samples-(chan_max+kk/bins))/(2*bwmul));
    end
    H = zeros(2*numel(samples)-1,1);
    H(1:numel(samples)) = sqrt(sum(abs(FILTS.^2),1));
    H(numel(samples)+1:end) = H(numel(samples)-1:-1:1); 

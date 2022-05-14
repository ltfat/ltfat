function H = comp_zerofilt(wintype,fs,chan_min,freqtoscale,scaletofreq,bwmul,bins,Ls)
%COMP_ZEROFILT low-pass filter for warped filter banks

%   see blfilter for an extended version
    kk = chan_min;
    while scaletofreq(kk+bwmul) > fs/Ls 
      kk = kk-1/bins;
    end
    Minfilt = kk;
    
    Maxpos = floor(Ls/fs*scaletofreq(chan_min-1/bins+bwmul));
    samples = freqtoscale((0:Maxpos)*fs/Ls);
    if samples(1) == -Inf
        samples(1) = samples(2);
    end
    
    FILTS = zeros(round(bins*(chan_min-Minfilt)),numel(samples));
    for kk = 1:size(FILTS,1)
       FILTS(kk,:) = firwin(wintype,(samples-(chan_min-kk/bins))/(2*bwmul));
    end
    H = zeros(2*numel(samples)-1,1);
    H(numel(samples):end) = sqrt(sum(abs(FILTS.^2),1));
    H(1:numel(samples)-1) = H(end:-1:numel(samples)+1); 

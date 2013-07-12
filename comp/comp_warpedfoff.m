function foff=comp_warpedfoff(fc,bw,fs,L,scaletofreq)
%COMP_WARPEDFOFF  foff for warped filters

foff=floor(scaletofreq(fc-.5*bw)/fs*L)+1;

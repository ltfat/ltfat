function foff=comp_warpedfoff(fc,bw,fs,L,freqtoscale,scaletofreq,do_symmetric)
%COMP_WARPEDFOFF  foff for warped filters

fcwasnegative = fc < 0;

if fcwasnegative && do_symmetric
   fc = -fc;
   fcscale = freqtoscale(fc);
   foff = -floor(scaletofreq(fcscale+.5*bw)/fs*L)+1;
else
   fcscale = freqtoscale(fc);
   foff = floor(scaletofreq(fcscale-.5*bw)/fs*L)+1;
end

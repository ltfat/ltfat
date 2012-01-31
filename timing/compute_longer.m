%COMPUTE_LONGER  Vary the length of the transform
%
%  This script computes the running time for longer and longer
%  transforms. Use the script plot_longer to visualize the result.
%
%  All other parameters except L remain fixed. The window length for the
%  filter bank algorithm is also kept fixed.
%
a=40;
M=60; 
W=4;
nrep=20;
gl=2400;
bl=24000;

% Get test sizes. Use only test sizes from nextfastfft, as the others
% cause to large a variability.
[dummy,table]=nextfastfft(1);

system('rm longer_fb_real.log');
system('rm longer_fac_real.log');
system('rm longer_ola_real.log');


testrange=M*table(find(table<=1500));
idx=round(testrange/gl)==testrange/gl;
testrange=testrange(idx);
testrange

for ii=1:length(testrange)
  L=testrange(ii);
  
  s=sprintf('./time_dgtreal_fb %i %i %i %i %i %i >> longer_fb_real.log\n',a,M,L,W, ...
            gl,nrep);
  disp(s);
  system(s);

  s=sprintf('./time_dgtreal_fac %i %i %i %i %i >> longer_fac_real.log\n',a,M,L,W,nrep);  
  disp(s);
  system(s);

  % Extend L to multiple of bl for the OLA algorithm.
  Lola = ceil(L/bl)*bl;

  s=sprintf('./time_dgtreal_ola %i %i %i %i %i %i %i >> longer_ola_real.log\n',a,M,Lola,W,gl,bl,nrep);  
  disp(s);
  system(s);
  
end;





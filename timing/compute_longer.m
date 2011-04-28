%COMPUTE_LONGER  Vary the length of the transform
%
%  This script computes the running time for longer and longer
%  transforms. Use the script plot_longer to visualize the result.
%
%  All other parameters except L remain fixed. The window length for the
%  filter bank algorithm is also kept fixed.
%
a=50;
M=200; 
W=4;
nrep=20;
gl=10*M;

% Get test sizes. Use only test sizes from nextfastfft, as the others
% cause to large a variability.
[dummy,table]=nextfastfft(1);

system('rm longer_fb1.log');
system('rm longer_fac1.log');
for L=gl:M:M*table(find(table<=500));
  s=sprintf('./time_dgt_fb %i %i %i %i %i %i >> longer_fb.log\n',a,M,L,W, ...
            gl,nrep);
  disp(s);
  system(s);

  s=sprintf('./time_dgt_fac %i %i %i %i %i >> longer_fac.log\n',a,M,L,W,nrep);  
  disp(s);
  system(s);
end;


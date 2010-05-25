%COMPUTE_LONGER  Vary the length of the transform
%
%  This script computes the running time for longer and longer
%  transforms. Use the script plot_longer to visualize the result.
%
%  All other parameters except L remain fixed. The window length for the
%  filter bank algorithm is also kept fixed.
%
a=320;
M=640; 
W=4;
nrep=20;
gl=10*M;

system('rm longer_fb.log');
system('rm longer_fac.log');
for L=gl:M:100*M
  s=sprintf('./time_dgt_fb %i %i %i %i %i %i >> longer_fb.log\n',a,M,L,W, ...
            gl,nrep);
  disp(s);
  system(s);

  s=sprintf('./time_dgt_fac %i %i %i %i %i >> longer_fac.log\n',a,M,L,W,nrep);  
  disp(s);
  system(s);
end;

%COMPUTE_LONGER  Vary the length of the transform
%
%  This script computes the running time for longer and longer
%  transforms. Use the script plot_longer to visualize the result.
%
%  All other parameters except L remain fixed. The window length for the
%  filter bank algorithm is also kept fixed.
%
a_init=16;
W=4;
nrep=20;

for a=a_init:a_init:50*a_init
  M=2*a;
  gl=10*a;
  for L=gl:M:50*M
    s=sprintf('./time_dgt_fb %i %i %i %i %i %i >> big_fb.log\n',a,M,L,W, ...
              gl,nrep);
    disp(s);
    system(s);
    
    s=sprintf('./time_dgt_fac %i %i %i %i %i >> big_fac.log\n',a,M,L,W,nrep);  
    disp(s);
    system(s);
  end;
end;



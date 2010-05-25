a=32;
M=64; 
L=a*M*100; 
W=4;
nrep=4;


system('rm crossover.log');
for gl=M:M:20*M
  s=sprintf('./time_dgt_fb %i %i %i %i %i %i >> crossover.log\n',a,M,L,W,gl,nrep);
  
  disp(s);
  system(s);
end;

s=sprintf('./time_dgt_fac %i %i %i %i %i > crossover.ref\n',a,M,L,W,nrep);

disp(s);
system(s);
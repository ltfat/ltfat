if 1
  % This is the setup used in the paper
  a=30;
  M=60; 
  L=a*M; 
  W=4;
  nrep=20;
else
  a=16;
  M=64; 
  L=a*M; 
  W=1;
  nrep=4;  
  
end;

system('rm crossover.log');
%for gl=M:M:20*M
for gl=M:M:16*M
  s=sprintf('./time_dgt_fb %i %i %i %i %i %i >> crossover.log\n',a,M,L,W,gl,nrep);
  
  disp(s);
  system(s);
end;

s=sprintf('./time_dgt_fac %i %i %i %i %i > crossover.ref\n',a,M,L,W,nrep);

disp(s);
system(s);